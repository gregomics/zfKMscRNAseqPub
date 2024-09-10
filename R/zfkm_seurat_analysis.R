library(tidyverse)
library(Seurat)
library(sctransform)
library(gridExtra)
library(scCATCH)
library(readxl)
library(pals)
library(ggalluvial)

CellRangeRdir <- "/CellRangeR/output/dir"


  
### Creating Seurat object

#Seurat (v 4.2.0) is used to analyse the scRNA-seq.
#First we need to create a Seurat object per sample and then merge them together.

# but first need to retrieve each individual seurat object.
# not that each seurat object is from an individual biological rep
getSeuratObj <- function(sampleID, dir, genotype) {
  datadir <- paste0(dir, "/", sampleID, "/")
  seurat.data <- Read10X(data.dir = datadir)
  seurat.obj <- CreateSeuratObject(counts = seurat.data,
                                   project = sampleID,
                                   min.cells = 3,
                                   min.features = 200)
  seurat.obj[["genotype"]] <- genotype
  return(seurat.obj)
}
wt_1.obj <- getSeuratObj("WT-1", CellRangeRdir, "WT")
wt_2.obj <- getSeuratObj("WT-2", CellRangeRdir, "WT")
rad21_1.obj <- getSeuratObj("RAD21-1", CellRangeRdir, "RAD21")
rad21_2.obj <- getSeuratObj("RAD21-2", CellRangeRdir, "RAD21")



all <- merge(wt_1.obj, y = c(wt_2.obj, rad21_1.obj, rad21_2.obj),
             add.cell.ids = c("wtr1", "wtr2", "rad21r1", "rad21r2"),
             project = "all")

### Exclusion of doublets:

#Cell doublets were estimated separately using scDblFinder (doi:10.12688/f1000research.73600.2).

#The script used to generate the scDblFinder output file can be found there: R/zfKM_scRNAseq_doublet_estimation.R

#Here we are parsing the scDblFinder output (1 file/sample). Then we will add this information in the meta.data slot from Seurat object.

sampleIDs <- c("WT-1", "WT-2", "RAD21-1", "RAD21-2")

scDblFinder.df <- data.frame()

for (sampleID in sampleIDs) {
  file <- paste0(CellRangeRdir, sampleID, "_scDblFinder.tsv")
  if (file.exists(file)) {
    summary <- read.table(file, sep = " ", header = T, stringsAsFactors = F)
    summary$sampleID <- sampleID 
    summary$cellBarcode <- paste0(gsub("-", "r",tolower(sampleID) ) , "_",
                                  summary$Barcode)
    scDblFinder.df <- rbind(scDblFinder.df, summary)
  } else {
    stop(paste0("Missing scDblFinder file: ", file))
  }

}

# plotting the doublets estimation
scDblFinder.df %>% 
  group_by(scDblFinder.class, sampleID) %>% 
  summarise(count = n()) %>% 
  ggplot(., aes(scDblFinder.class, count, fill = scDblFinder.class)) + 
  geom_bar(stat = "identity") + 
  facet_grid(~sampleID) + 
  theme_bw() + 
  ggtitle("Doublet estimation using scDblFinder")

rownames(scDblFinder.df) <- scDblFinder.df$cellBarcode

# now adding the doubletstatus column in the seurat object:

all@meta.data$doubletstatus <- scDblFinder.df[rownames(all@meta.data),]$scDblFinder.class

#all@meta.data %>% ggplot(., aes(nCount_RNA, nFeature_RNA, col = doubletstatus)) + geom_point(alpha = 0.5) + facet_grid(~orig.ident) + theme_bw() + ggtitle("Doubet estimation using scDblFinder")

### QC from Seurat object

#QC data are stored into meta.data object.

# we can add more data in QC using Seurat functions:
all[["percent.mt"]] <- PercentageFeatureSet(all, pattern = "^mt-")

### Filtering using scDblFinder and percentage mitochondrial transcript < 5%

all<- subset(all, subset = doubletstatus == "singlet" & percent.mt < 5)

#Data will then be normalised by default Seurat uses a lognormal approach and it further use a scaling factor of 1000.
#Once the data is normalised, then variable genes will be identified using the vst method. 


## creating the list of Seurat object. split by genotype.
geno.list <- SplitObject(all, split.by = "genotype")

# for all seurat object in the list: Norm + FindVariable feature
geno.list <- lapply(X = geno.list, FUN = function(x) {
  x <- NormalizeData(x)  
  x <- FindVariableFeatures(x, selection.method = "vst", 
                            nfeatures = 2000)
})


### Now integrating all together
## meaning Seurat will identify "anchors" correspondences (genes with same profile/cell) between conditions. 
## use these anchors to harmonise the 2 dataset
# more info there: https://doi.org/10.1016/j.cell.2019.05.031
geno.anchors <- FindIntegrationAnchors(object.list = geno.list, 
                                       dims = 1:20)
## This step is time consuming. ~1hr with 12Gb RAM

geno.combined <- IntegrateData(anchorset = geno.anchors, dims = 1:20)


#After this step, the Seurat object (geno.combined) will have 2 assays: integrated and RNA. Integrated assay should be used for the PCA, UMAP and cluster formation. RNA assay should be used for differential gene expression analysis purpose.


DefaultAssay(geno.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
geno.combined <- ScaleData(geno.combined, verbose = FALSE)
geno.combined <- RunPCA(geno.combined, npcs = 30,
                        verbose = FALSE)
# UMAP and Clustering
geno.combined <- RunUMAP(geno.combined, reduction = "pca",
                         dims = 1:20)
geno.combined <- FindNeighbors(geno.combined,
                               reduction = "pca", dims = 1:20)

## here we set the resolution to 0.5
geno.combined <- FindClusters(geno.combined, resolution = 0.5)

DimPlot(geno.combined, reduction = "umap",
        split.by = "genotype", label = TRUE) +  ggtitle("UMAP. Res:0.5") + theme(text = element_text(size=8))

geno.clustIdents <- Idents(geno.combined)


### Identitifying conserved markers between condition by cluster

#We need to loop to the identified clusters and check what are the differentially expressed genes (markers) from that cluster compared to the other clusters that are conserved between condition.

# extracting clusters:
clusterIdents <- levels(geno.clustIdents)

# FindConservedMarkers will look at differential marker from condition 1 in clusterA against all the other clusters.

# looping through clusters. step is time consuming: ~1hr
conservedMarks.df <- data.frame()
for (cl in clusterIdents) {
  message(paste0("processing cluster: ", cl, "<br>" ))
  comMarks.file <-paste0(CellRangeRdir, "conservedMarkers_", cl, "_res0.5.xls")
  if (file.exists(comMarks.file)) {
    message(paste0(comMarks.file, " is already there!\n"))
  }
  else {
    comMark <- FindConservedMarkers(geno.combined, ident.1 = cl, 
                                    grouping.var = "genotype", verbose = T)
    comMark %>% as.data.frame() %>% rownames_to_column() %>%  write.table(., comMarks.file, sep = "\t", quote = F)
  }
  comMarks.data <- read.table(comMarks.file,
                                sep = "\t",
                                stringsAsFactors = F,
                                header = T)
  comMarks.data$cluster <- cl
  conservedMarks.df <- rbind(conservedMarks.df, comMarks.data)
}


### Differentially Expressed Gene between Rad21a+/- vs WT by cluster

#To retrieve genes that are differentially expressed between Rad21a+/- vs WT, we used the DESeq2 method in the Seurat::FindMarkers function.
#Prior to that, it is important to set the DefaultAssay as the RNA slot.


Idents(geno.combined) <- "celltype.geno"

RNA_dges.df <- data.frame()
DefaultAssay(geno.combined) <- "RNA"
for (cl in levels(geno.combined@meta.data$celltype) ) {
  cat(paste0("Processing cluster: ", cl, "\n"))
  # processing RAD21 v WT
  grp1 <- paste0(cl, "_RAD21")
  grp2 <- paste0(cl, "_WT")
  clustdge_rad21_v_WT <- FindMarkers(geno.combined, 
                                     ident.1 = grp1, ident.2 = grp2,
                                     verbose = T,
                                     test.use =  "DESeq2",
                                     slot = "counts") %>% 
    as.data.frame() %>% 
    rownames_to_column()
    
  if ( nrow(clustdge_rad21_v_WT) > 0 ) {
    clustdge_rad21_v_WT$cluster <- cl
    clustdge_rad21_v_WT$comparison <- "Rad21_v_WT"
    RNA_dges.df <- rbind(RNA_dges.df, clustdge_rad21_v_WT )
  }
    
  
}
# setting back to the default assay
DefaultAssay(geno.combined) <- "integrated"

RNA_dges.file <-paste0(CellRangeRdir,
                       "DESeq_DGEs_between_conditions_res0.5_all_new.xls")
RNA_dges.df %>% write.table(., RNA_dges.file, sep = "\t", quote = F)


### Automatic annotation using scCATCH

#scCATCH is an automatic tool that assign cellular identity. It is reference based meaning that it uses a repository of curated data. This repository is named cellMatch
#To create a scCATCH object we need a sparse matrix and a vector of clusterID. scCATCH is normally used for Human, so we used the orthologuous genes retrieve from martview (Ensembl).

zf2hsOrthologFile <- paste0(CellRangeRdir,
                            "/martview_zf_v_hs_ortolog.txt")
if (file.exists(zf2hsOrthologFile)) {
  zf2hsOrtholog <- read.table(zf2hsOrthologFile,
                              sep = "\t",
                              header = T,
                              stringsAsFactors = F)  
} else {
  cat(paste0("missing ortholog file: ", zf2hsOrthologFile))
}

# extracting the sparse matrix
scCATCH.obj <- createscCATCH(data = geno.combined@assays$RNA@data, 
                             cluster = as.character(geno.combined$celltype))
### need to construct custom cellmarker to be use to assign cell identity.
hs.cellmatch <- cellmatch %>% 
  dplyr::filter(species == "Human") %>% 
  dplyr::filter(   cancer == "Normal")
zf.cellmatch <- hs.cellmatch %>% 
  left_join(., zf2hsOrtholog %>% dplyr::select(c(HumanGeneName, zfGeneName)), 
            by = c("gene" = "HumanGeneName")) %>%  
  dplyr::select(-c(gene, species)) %>% 
  dplyr::filter(! is.na(zfGeneName )) %>% 
  rename(gene = zfGeneName) %>% 
  mutate(species = "Zebrafish")

scCATCH.obj <- findmarkergene(object = scCATCH.obj,
                              marker = zf.cellmatch,
                              if_use_custom_marker = TRUE)

scCATCH.obj <- findcelltype(scCATCH.obj)
scCATCH.obj@celltype 


### CellIdentity using data from publication (JEM_20170976_TableS3)

#Table S3 from the publication contains cell identity annotation.
#We need to convert that format to the one fitting with the scCATCH tool.
#The main issue is that each cell type data is contained into an excel sheet.
#So first need to extract this info.   

### structure of the cellmatch file:
# colnames(zf.cellmatch)
# [1] "tissue"    "cancer"    "condition" "subtype1"  "subtype2"  "subtype3"  "celltype" 
# [8] "resource"  "pmid"      "gene"      "species"

## looking at the zf.cellmatch overlapping genes being in several different cell type is allowed.

jem_path <- paste0(CellRangeRdir, "JEM_20170976_TableS3.xlsx")
## Using read_excel from library(readxl) we will be able to extract it:
jem_data <- jem_path %>%  excel_sheets() %>% set_names() %>%  map(read_excel, path = jem_path)
#jem_data will be a list with the name being the sheet name. 

jem.df <- data.frame()
for (identity in names(jem_data)) {
  identity.df <- jem_data[[identity]]
  cat(paste0("Processing ", identity, " nrow: ", nrow(identity.df), "\n"))
  tmp.df <-data.frame(gene = identity.df$genename,
                      condition = "Normal",
                      celltype = identity,
                      cancer = "Normal",
                      tissue = "Kidney",
                      subtype1 = NA,
                      subtype2 = NA,
                      subtype3 = NA,
                      resource = "Experiment",
                      pmid = "20170976",
                      species = "Zebrafish")
  jem.df <- rbind(jem.df, tmp.df)
}

out.file <- paste0(CellRangeRdir, "JEM_20170976_TableS3_cellmatch.tsv")
jem.df %>% write.table(.,out.file, quote = F, sep = "\t")


### Manual annotation

#After manual inspection and curation of cell type, we provide a cell annotation.

annotation_file <- paste0(CellRangeRdir, "final_cluster_annotation.tsv")
manual_annotation <- read.table(annotation_file,
                                sep = "\t",
                                header = T,
                                stringsAsFactors = F)
cluster_annotation <- geno.combined@meta.data %>% 
  mutate(cluster = as.integer(as.character(seurat_clusters))) %>%  
  left_join(., manual_annotation, by = c("cluster" = "cluster_no")) %>% 
  pull(annotation)
# adding to the metadata 

geno.combined@meta.data$cluster_annotation <- cluster_annotation

saveRDS(geno.combined, "geno.combined.rds")

