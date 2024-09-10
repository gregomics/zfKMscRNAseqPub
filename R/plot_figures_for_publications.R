## to plot scRNA-seq figures for the publications, we used a preprocessed RDS object geno.combined
# from there we load the data and uses either Seurat, monocle3 or ggplot2 plotting function.
library(Seurat)
library(tidyverse)
library(monocle3)
library(SeuratWrappers) # use to convert seurat to sce


seuratdir <- "/path_to_where_is_stored_rds_file"
filename <- "geno.combined.rds"
rds.file <- paste0(seuratdir, filename)
geno.combined <- readRDS(rds.file)


### Figure1A: This is the plotMA from the bulk RNA-seq.

### Figure 2A: 

## plot high resolution UMAP for manuscript. We need to change heatmap

### high resolution UMAP with a different color schema:
# to plot the UMAP plot of each genotype we need to split the geno.combined object by genotype:

seuratobj.geno.list <-  SplitObject(geno.combined, split.by = "genotype")
rad21seurat.obj <- seuratobj.geno.list[['RAD21']]
wtseurat.obj <- seuratobj.geno.list[['WT']]

### WT UMAP plot:

tiff.file <- paste0(seuratdir, "UMAP_wt_only_190324.tiff")
tiff(tiff.file, res = 600, height = 20, width = 20, units = "cm")
DimPlot(wtseurat.obj, reduction = "umap",
        label = F, cols = c('0' = "#E16272", '1' = "#CF2B57", '2' = "#6DC7B4", '3' = "#4D739F", '4'="#225FAC", '5'="#46B449", '6' = "#8FB2DA", '7' = "#F7941D", '8' = "#CD81B7", '9' = "#92278F", '10' = "#AC98C9", '11' = "#9B9B9B", '12' = "#F7941D", '13' = "#9B9B9B", '14' = "#84C441", '15' = "#3986C7", '16' =  "#FFF200", '17' = "#F7941D", '18' = "#9B9B9B", '19' = "#9B9B9B", '20' = "#F7941D", '21'="#9B9B9B")) +  ggtitle("UMAP. WT genotype only") + theme(text = element_text(size=8))
dev.off()

### Rad21 ko

tiff.file <- paste0(seuratdir, "UMAP_rad21_only_190324.tiff")
tiff(tiff.file, res = 600, height = 20, width = 20, units = "cm")
DimPlot(rad21seurat.obj, reduction = "umap",
        label = F, cols = c('0' = "#E16272", '1' = "#CF2B57", '2' = "#6DC7B4", '3' = "#4D739F", '4'="#225FAC", '5'="#46B449", '6' = "#8FB2DA", '7' = "#F7941D", '8' = "#CD81B7", '9' = "#92278F", '10' = "#AC98C9", '11' = "#9B9B9B", '12' = "#F7941D", '13' = "#9B9B9B", '14' = "#84C441", '15' = "#3986C7", '16' =  "#FFF200", '17' = "#F7941D", '18' = "#9B9B9B", '19' = "#9B9B9B", '20' = "#F7941D", '21'="#9B9B9B")) +  ggtitle("UMAP. RAD21KO genotype only") + theme(text = element_text(size=8))
dev.off()


### Figure 2B: Stacked bar of % cells population

clusters_colors <- c('0' = "#E16272", '1' = "#CF2B57", '2' = "#6DC7B4", '3' = "#4D739F", '4'="#225FAC", '5'="#46B449", '6' = "#8FB2DA", '7' = "#F7941D", '8' = "#CD81B7", '9' = "#92278F", '10' = "#AC98C9", '11' = "#9B9B9B", '12' = "#F7941D", '13' = "#9B9B9B", '14' = "#84C441", '15' = "#3986C7", '16' =  "#FFF200", '17' = "#F7941D", '18' = "#9B9B9B", '19' = "#9B9B9B", '20' = "#F7941D", '21'="#9B9B9B")

geno.combined@meta.data  %>% dplyr::filter(genotype %in% c("RAD21", "WT")) %>% group_by(orig.ident, seurat_clusters) %>% summarise(cluster_count = n()) %>% dplyr::left_join(., geno.combined@meta.data  %>% dplyr::filter(genotype %in% c("RAD21", "WT")) %>% group_by(orig.ident) %>% summarise(total_count = n()) , by = c("orig.ident" = "orig.ident")) %>% dplyr::mutate(clust_perc = (cluster_count / total_count) * 100) %>% ggplot(., aes(orig.ident, clust_perc, fill = seurat_clusters)) + geom_bar(stat = "identity") + theme_bw() + scale_fill_manual(values = clusters_colors)


### Figure 2C. 
# This figure was done in PRISM. This is the R equivalent.

geno.combined@meta.data  %>% 
  dplyr::filter(genotype %in% c( "WT")) %>% 
  group_by(genotype, seurat_clusters) %>% 
  summarise(wt_cluster_count = n()) %>% 
  dplyr::left_join(., geno.combined@meta.data  %>% 
                     dplyr::filter(genotype %in% c( "RAD21")) %>% group_by(genotype, seurat_clusters) %>% 
                     summarise(rad21_cluster_count = n()), by = c("seurat_clusters" = "seurat_clusters")) %>% 
  dplyr::mutate(fc = 1 - (rad21_cluster_count / wt_cluster_count)) %>% 
  ggplot(., aes(fc, seurat_clusters)) + geom_point() + theme_bw()

### Figure 2D:
# this is a composite of 3 expression profiles.

genes <- c( "pcna", "myca", "rad21a")
# first we need to extract nonkidney cells, the iskidney flag was manually curated
# based on the kidney gene markers.

nonkidney.obj <- subset(geno.combined, iskidney == "N")

nonkidney.obj <- subset(nonkidney.obj, genotype %in% c("WT", "RAD21"))

nonkidney.norm.exp <- nonkidney.obj[["RNA"]]@data %>% as.data.frame()

for (gene in genes) {
  nonkidney.norm.exp.gene.df <- nonkidney.norm.exp %>% rownames_to_column("genename") %>% dplyr::filter(genename == gene ) %>% pivot_longer(!genename, names_to =  "BC") %>% dplyr::left_join(., nonkidney.obj@meta.data %>%  rownames_to_column() %>% dplyr::select(c(rowname, seurat_clusters, genotype)) , by = c("BC" = "rowname")) %>% mutate(cluster = factor(seurat_clusters, levels = clusters_order))
  title <- paste0(gene, " expression profile.")
  q <- nonkidney.norm.exp.gene.df %>% 
    ggplot(., aes(cluster, value, fill = genotype)) + 
    geom_violin(scale = "width", trim = 0.2 )  + 
    theme_bw() + 
    geom_pointrange(mapping = aes(x = cluster, y = value), 
                    stat = "summary", fun.min = function(x) {quantile(x,0.25)}, 
                    fun.max = function(x) {quantile(x,0.75)}, fun="median", show.legend = F, 
                    position = position_dodge(.800), size = 0.1, col = "darkblue") + 
    ggtitle(title)
  pdf.file <- paste0(seuratdir, "nonkidney_", gene, "_violin_pointrange.pdf")
  pdf(pdf.file)
  print(q)
  dev.off()
  
}

### plotting actin genes: actb1, actb2 along erythroids cluster

eryt_clusters    <- c("8","9","0","1")
erythroids  <- subset(geno.combined, 
                      subset = seurat_clusters %in% eryt_clusters & genotype %in% c("WT", "RAD21"))
erythroids.norm.exp <- erythroids[["RNA"]]@data %>% as.data.frame()
actin_genes <- c("actb1", "actb2")

wt_eryt_order <- c("8", "9", "0", "1")

for (gene in actin_genes) {
  erythroids.norm.exp.gene.df <- erythroids.norm.exp %>% 
    rownames_to_column("genename") %>% 
    dplyr::filter(genename == gene ) %>% 
    pivot_longer(!genename, names_to =  "BC") %>% 
    dplyr::left_join(.,  erythroids@meta.data %>%  
                       rownames_to_column() %>% dplyr::select(c(rowname, seurat_clusters, genotype)) , 
                     by = c("BC" = "rowname")) %>% 
    mutate(cluster = factor(seurat_clusters, levels = wt_eryt_order))
  
  title <- paste0(gene, " expression profile in Erythrocytes.")
  q <- erythroids.norm.exp.gene.df %>% 
    ggplot(., aes(cluster, value, fill = genotype)) + 
    geom_violin(scale = "width", trim = 0.2 )  + 
    theme_bw() + 
    geom_pointrange(mapping = aes(x = cluster, y = value), 
                    stat = "summary", fun.min = function(x) {quantile(x,0.25)}, 
                    fun.max = function(x) {quantile(x,0.75)}, fun="median", 
                    show.legend = F, position = position_dodge(.800), size = 0.1, col = "darkblue") + 
    ggtitle(title)
  
  pdf.file <- paste0(seuratdir, "eryt_", gene, "_violin_pointrange.pdf")
  pdf(pdf.file)
  print(q)
  dev.off()
  
}



### Figure 2E. Done using PRISM. We could generate similar plot using R and the DESeq_DGEs_between_conditions_res0.5_all.xlsx
# which is the output genotype gene expression comparison at each cluster.
# note this is one of the file generated by Seurat processing/Analysis file.

RNA_dges.file <-paste0(seuratdir,
                       "DESeq_DGEs_between_conditions_res0.5_all.xls")
RNA_dges.df <- read.table(RNA_dges.file, sep = "\t", header = T,
                          stringsAsFactors = F)

## Figure 2F. MPP Gene Ontology pathway.
## this figure is build using RNA_dges.df dataframe containing differentially expressed genes using DESeq2.
## it will use DOSE, ClusterProfiler and OrgDB for Symbol to ID conversion.

library(DOSE)
library(ClusterProfiler)
library(OrgDB)

clust <- 10
cutoff <- 0.5
mpp_sign_dges <- RNA_dges.df %>% 
  dplyr::filter(comparison == "Rad21_v_WT", p_val_adj <0.05, cluster == clust ) %>%  
  mutate(reg = ifelse(avg_log2FC >0 , "up", "down"))

# getting ENTREZID
gene.df <- bitr(mpp_sign_dges$rowname,
                fromType = "SYMBOL",
                toType = "ENTREZID",
                OrgDb = OrgDB)
# extracting MPP differentially expressed genes:
mpp_sign_dges <- mpp_sign_dges %>% 
  dplyr::left_join(., gene.df, by=c("rowname" = "SYMBOL"))

# comparing BP GO enrichment in function of the regulation (up/down)
mpp_sign_dges.ck <- compareCluster(ENTREZID~reg, data=mpp_sign_dges,
                                   fun="enrichGO", OrgDb = OrgDB, 
                                   ont = "BP",
                                   readable = T)
# to avoid replication of BP gene ontology having nearly similar content
# we simplified it using the clusterProfiler::simplify at the highest level

mpp_sign_dges.ck.simp <- clusterProfiler::simplify(mpp_sign_dges.ck,
                                                   cutoff=cutoff,
                                                   by="p.adjust",
                                                   select_fun = min)

## plotting the results:
mpp_sign_dges.ck.simp %>% 
  as.data.frame() %>%   
  mutate(regF = factor(reg, levels = c("up", "down"))) %>%    
  ggplot(., aes(regF, Description, col = -log10(p.adjust), size = Count)) + 
  geom_point() + 
  scale_color_continuous(type = "viridis") + 
  theme_bw()  + 
  theme(axis.text.y.left = element_text(size = 6)) + 
  ggtitle("ClusPro. Simp 0.5. Most generic terms")

### end of figure 2F.

### figure 3A: UMAP and monocle3 trajectory analysis:
## trajectory analysis will be performed on erythrocyte only using MPP (seurat_clusters 10) as the root.

library(SeuratWrappers)
library(monocle3)
library(pheatmap)
library(RColorBrewer)
library(rstatix)
library(ggridges)

working_dir <- ""
rds_file <- paste0(working_dir, "geno.combined.rds")
geno.combined <-readRDS(rds_file)
clusters_cols = c('0' = "#E16272", '1' = "#CF2B57", '2' = "#6DC7B4", '3' = "#4D739F", '4'="#225FAC", '5'="#46B449", '6' = "#8FB2DA", '7' = "#F7941D", '8' = "#CD81B7", '9' = "#92278F", '10' = "#AC98C9", '11' = "#9B9B9B", '12' = "#F7941D", '13' = "#9B9B9B", '14' = "#84C441", '15' = "#3986C7", '16' =  "#FFF200", '17' = "#F7941D", '18' = "#9B9B9B", '19' = "#9B9B9B", '20' = "#F7941D", '21'="#9B9B9B")


### Erythrocyte only using MPP as root. 

eryt_clusters    <- c(10,8,9,0,1)
erythroids  <- subset(geno.combined, 
                      subset = seurat_clusters %in% eryt_clusters & genotype %in% c("WT", "RAD21"))

# conversion from Seurat object to cell_data_set object
erythroids.cds <- as.cell_data_set(erythroids)

## creating the list of partitions. 
erythroids.partitions <- rep(1, length(erythroids.cds@colData@rownames))

names(erythroids.partitions) <- erythroids.cds@colData@rownames
erythroids.partitions <- factor(erythroids.partitions)
erythroids.cds@clusters$UMAP$partitions <- erythroids.partitions

### Now assigning the clusters info:
erythroids.cds@clusters$UMAP$clusters <- erythroids@active.ident

## Need to add the UMAP embedding coordinates from Seurat
# the dimensional reductions are stored in the cds object there:
#reducedDimNames(2): PCA UMAP
# normally this is already stored
erythroids.cds@int_colData@listData$reducedDims$UMAP <- erythroids@reductions$umap@cell.embeddings

### Adding gene names to the object as it seems that seuratWrappers does not propagate it.

erythroids.cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(geno.combined[["RNA"]])
# this fixes the issue of plot_cells with a list of gene names

#using the utility as.cell_data_set create by default the "Size_Factor" column in the colData dataframe. The issue is that column does not really contains a Size_Factor but the nCount.
erythroids.cds <- estimate_size_factors(erythroids.cds)
### According to monocle3 code:
#cell_total <- Matrix::colSums(counts)
#sfs <- cell_total / exp(mean(log(cell_total)))

## We cannot use the partition because that was not calculated here. Instead we 
erythroids.cds <- learn_graph(erythroids.cds,
                              use_partition = F)


### reordering the cells based on the pseudotime.
### It needs a preliminary knowledge such as the roots.
### here we can start by setting the cluster 10 (MPP) as a root.

cluster10cells <- colData(erythroids.cds) %>% 
  as.data.frame() %>% 
  rownames_to_column("cellids") %>% 
  dplyr::filter(celltype == 10) %>% 
  pull(cellids) 

erythroids.cds <- order_cells(erythroids.cds,
                              reduction_method = "UMAP",
                              root_cells = cluster10cells)
## Note that the pseudotime will be stored there for each cell:
#erythroids.cds@principal_graph_aux@listData$UMAP$pseudotime

q_eryt_pseudo <- plot_cells(erythroids.cds,
                            color_cells_by = "pseudotime", 
                            label_branch_points = F, 
                            label_leaves = F,
                            label_roots = F,
                            trajectory_graph_color = "green") + theme(legend.position = 'top')


eryt_colors <- clusters_cols[names(clusters_cols) %in% eryt_clusters] 

q_eryt_cluster <- plot_cells(erythroids.cds,
                             color_cells_by = "seurat_clusters", 
                             label_branch_points = F, 
                             label_leaves = F,
                             label_roots = F,
                             label_cell_groups = F,
                             label_groups_by_cluster = F,
                             trajectory_graph_color = "green") + scale_color_manual(values = eryt_colors) + theme(legend.position = 'top')

pdf(paste0(working_dir, "erythrocyte_monocle3_2008.pdf"))
grid.arrange(q_eryt_pseudo, q_eryt_cluster)
dev.off()

### end of figure 3A

### figure 3B: Erythrocyte violin plots per cluster and genotype:

erythroids.cds$monocle_pseudotime <- pseudotime(erythroids.cds)
erythroids.df <- colData(erythroids.cds) %>%  
  as.data.frame()
erythroids.df %>% 
  dplyr::filter(genotype == "WT") %>% 
  ggplot2::ggplot(., aes(monocle_pseudotime,reorder( seurat_clusters, monocle_pseudotime, median), 
                         fill = genotype)) + 
  geom_jitter(aes(monocle_pseudotime,reorder( seurat_clusters, monocle_pseudotime, median), 
                  col = genotype),  alpha=0.5, size =0.5, position = position_jitterdodge()) + 
  geom_violin(scale="width", alpha = 0.5) + 
  geom_pointrange(mapping=aes(monocle_pseudotime,reorder( seurat_clusters, monocle_pseudotime, median)),
                  stat = "summary", fun.min = function(x) {quantile(x,0.25)}, 
                  fun.max = function(x) {quantile(x,0.75)}, 
                  fun="median", show.legend = F, size = 0.3, 
                  col = "darkblue", position=position_dodge(width=1))  + 
  theme_bw() + 
  ggtitle("WT Eryt processed together\nShowing both genotypes monocole2 pseudotime") + 
  ylab("cell type")

wt_eryt_order <- c("10", "8", "9", "0", "1")
rad21_eryt_order <- c("10", "9", "8", "1", "0")
# to be consistent with the other figure, it needs to be colored accordingly.
# ggplot2 uses hue_pal() available from scales package.
wt_rad21_hue <- hue_pal()(2)

wt_eryt_viol <- erythroids.df %>% 
  dplyr::filter(genotype == "WT") %>% 
  dplyr::mutate(cluster = factor(seurat_clusters, levels = wt_eryt_order)) %>% 
  ggplot2::ggplot(., aes(cluster, monocle_pseudotime, fill = genotype)) + 
  geom_jitter(aes(cluster, monocle_pseudotime, col = genotype),  
              alpha=0.5, size =0.5, position = position_jitterdodge()) + 
  geom_violin(scale="width", alpha = 0.5) + 
  geom_pointrange(mapping=aes(cluster, monocle_pseudotime), stat = "summary", 
                  fun.min = function(x) {quantile(x,0.25)}, 
                  fun.max = function(x) {quantile(x,0.75)}, 
                  fun="median", show.legend = F, size = 0.3, col = "darkblue", 
                  position=position_dodge(width=1))  + 
  theme_bw() + 
  ggtitle("WT Eryt processed together") + 
  ylab("cell type") + 
  scale_fill_manual(values = wt_rad21_hue[2]) + 
  scale_color_manual(values = wt_rad21_hue[2])

rad21_eryt_viol <- erythroids.df %>% 
  dplyr::filter(genotype == "RAD21") %>% 
  dplyr::mutate(cluster = factor(seurat_clusters, levels = rad21_eryt_order)) %>% 
  ggplot2::ggplot(., aes(cluster, monocle_pseudotime, fill = genotype)) + 
  geom_jitter(aes(cluster, monocle_pseudotime, col = genotype),  
              alpha=0.5, size =0.5, position = position_jitterdodge()) + 
  geom_violin(scale="width", alpha = 0.5) + 
  geom_pointrange(mapping=aes(cluster, monocle_pseudotime), stat = "summary", 
                  fun.min = function(x) {quantile(x,0.25)}, 
                  fun.max = function(x) {quantile(x,0.75)}, 
                  fun="median", show.legend = F, size = 0.3, 
                  col = "darkblue", position=position_dodge(width=1))  + 
  theme_bw() + 
  ggtitle("RAD21 Eryt processed together") + 
  ylab("cell type") + 
  scale_fill_manual(values = wt_rad21_hue[1]) + 
  scale_color_manual(values = wt_rad21_hue[1])

pdf(paste0(working_dir, "eryt_rad21_wt_pseudo_violin_2008.pdf"))
grid.arrange(wt_eryt_viol, rad21_eryt_viol)
dev.off()


### end of figure 3B

### figure figure 3C
# this figure was done using PRISM on DESeq2 DEG.

### end of figure 3C

### Figure 3D

genes <- c("gata1a", "alas2", "cahz", "hbaa1", "hbaa2",
           "hbba1", "gata1a", "znfl2a", "klf1",
           "hemgn", "si:ch1073-184j22.1")
clusters <- c(8,9,0,1 )

eryt.obj <- subset(geno.combined, 
                   subset = seurat_clusters %in% clusters & genotype %in% c("WT", "RAD21") )

# extracting the normalised count:

eryt_norm_exp <- GetAssayData(object = eryt.obj,
                              assay = "RNA",
                              slot = "data")

for (gene in genes ) {
  message(paste0("processing: ", gene))
  
  gene_norm_exp <- eryt_norm_exp[gene,] %>% 
    as.data.frame() %>% 
    dplyr::rename(value = ".") %>% 
    dplyr::mutate(gene = gene) %>% 
    rownames_to_column("BC")
  
  norm_exp.df <- eryt.obj@meta.data %>% 
    as.data.frame() %>%  
    rownames_to_column() %>% 
    dplyr::select(c(rowname, genotype, seurat_clusters)) %>%  
    dplyr::left_join(., gene_norm_exp , by = c("rowname" = "BC"))
  
  p <- norm_exp.df %>% 
    mutate(cluster = factor(seurat_clusters, levels =  clusters)) %>% 
    ggplot(., aes(cluster, value, fill= genotype, col = genotype))  + 
    theme_bw() + 
    geom_violin(size = 0.5, scale = "width")  + 
    facet_grid(gene~.) + 
    geom_pointrange(mapping = aes(x = cluster, y = value), stat = "summary", 
                    fun.min = function(x) {quantile(x,0.25)}, 
                    fun.max = function(x) {quantile(x,0.75)}, 
                    fun="median", show.legend = F, size = 0.1, col = "darkblue", 
                    position = position_dodge(width = 0.9)) + 
    theme(strip.text.y.right = element_text(angle = 0))
  
  plot.pdf <- paste0(outdir,"eryt_", gene, "_violinplot_1207.pdf")
  pdf(plot.pdf)
  print(p)
  dev.off()
}

### end of figure 3D

### Figure 3E: Metascape enrichment for Macrophages
### end of figure 3E

### Figure 3F: This was done with PRISM using DESeq2 DEG for macrophages
### end of figure 3F

### Figure 4A: Granulocyte UMAP and monocle3 trajectory analysis

gran_clusters <-  c(10, 15, 6, 3, 4 )
granulocytes <- subset(geno.combined, 
                       subset = celltype %in% gran_clusters & genotype %in% c("RAD21", "WT"))


granulocytes.cds <- as.cell_data_set(granulocytes)

granulocytes.partitions <- rep(1, length(granulocytes.cds@colData@rownames))

names(granulocytes.partitions) <-granulocytes.cds@colData@rownames
granulocytes.partitions <- factor(granulocytes.partitions)
granulocytes.cds@clusters$UMAP$partitions <- granulocytes.partitions

### Now assigning the clusters info:
granulocytes.cds@clusters$UMAP$clusters <- granulocytes@active.ident

## Need to add the UMAP embedding coordinates
# the dimensional reductions are stored in the cds object there:
#reducedDimNames(2): PCA UMAP
# normally this is already stored
granulocytes.cds@int_colData@listData$reducedDims$UMAP <- granulocytes@reductions$umap@cell.embeddings

## Adding gene names to the object as it seems that seuratWrappers does not propagate it.

granulocytes.cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(geno.combined[["RNA"]])

## We cannot use the partition since we set them as false.
granulocytes.cds <- learn_graph(granulocytes.cds,
                                use_partition = F)

### reordering the cells based on the pseudotime.
### It needs a preliminary knowledge such as the roots.
### here we can start by setting the cluster 1 as a root.

cluster10cells <- colData(granulocytes.cds) %>% as.data.frame() %>% rownames_to_column("cellids") %>% dplyr::filter(celltype == 10) %>% pull(cellids) 

granulocytes.cds <- order_cells(granulocytes.cds,
                                reduction_method = "UMAP",
                                root_cells = cluster10cells)


q_gran_pseudo <- plot_cells(granulocytes.cds,
                            color_cells_by = "pseudotime", 
                            label_branch_points = F, 
                            label_leaves = F,
                            label_roots = F,
                            trajectory_graph_color = "green") + theme(legend.position = 'top')


gran_colors    <- clusters_cols[names(clusters_cols) %in% gran_clusters] 

q_gran_cluster <- plot_cells(granulocytes.cds,
                             color_cells_by = "seurat_clusters", 
                             label_branch_points = F, 
                             label_leaves = F,
                             label_roots = F,
                             label_cell_groups = F,
                             label_groups_by_cluster = F,
                             trajectory_graph_color = "green") + 
  scale_color_manual(values = gran_colors) + 
  theme(legend.position = 'top')

pdf(paste0(working_dir, "granulocyte_monocle3_2008.pdf"))
grid.arrange(q_gran_pseudo, q_gran_cluster)
dev.off()

### end of figure 4A

### figure 4B: Granulocyte pseudotime violin plot.

granulocytes.cds$monocle_pseudotime <- pseudotime(granulocytes.cds)
granulocytes.df <- colData(granulocytes.cds) %>%  
  as.data.frame()
granulocytes.df %>% 
  dplyr::filter(genotype == "WT") %>% 
  ggplot2::ggplot(., aes(monocle_pseudotime,reorder( seurat_clusters, monocle_pseudotime, median), 
                         fill = genotype)) + 
  geom_jitter(aes(monocle_pseudotime,reorder( seurat_clusters, monocle_pseudotime, median), 
                  col = genotype),  alpha=0.5, size =0.5, position = position_jitterdodge()) + 
  geom_violin(scale="width", alpha = 0.5) + 
  geom_pointrange(mapping=aes(monocle_pseudotime,reorder( seurat_clusters, monocle_pseudotime, median)), 
                  stat = "summary", 
                  fun.min = function(x) {quantile(x,0.25)}, 
                  fun.max = function(x) {quantile(x,0.75)}, 
                  fun="median", show.legend = F, size = 0.3, 
                  col = "darkblue", position=position_dodge(width=1))  + 
  theme_bw() + 
  ggtitle("WT Eryt") + 
  ylab("cell type")

wt_gran_order <- c("10", "15", "3", "6", "4")
# here the order of granulocyte pseudotime median does not change between wt and rad21.
rad21_gran_order <- wt_gran_order
# to be consistent with the other figure, it needs to be colored accordingly.
# ggplot2 uses hue_pal() available from scales package.
wt_rad21_hue <- hue_pal()(2)

wt_gran_viol <- granulocytes.df %>% 
  dplyr::filter(genotype == "WT") %>% 
  dplyr::mutate(cluster = factor(seurat_clusters, levels = wt_gran_order)) %>% 
  ggplot2::ggplot(., aes(cluster, monocle_pseudotime, fill = genotype)) + 
  geom_jitter(aes(cluster, monocle_pseudotime, col = genotype),  alpha=0.5, size =0.5, 
              position = position_jitterdodge()) + 
  geom_violin(scale="width", alpha = 0.5) + 
  geom_pointrange(mapping=aes(cluster, monocle_pseudotime), stat = "summary", 
                  fun.min = function(x) {quantile(x,0.25)}, 
                  fun.max = function(x) {quantile(x,0.75)}, 
                  fun="median", show.legend = F, size = 0.3, 
                  col = "darkblue", position=position_dodge(width=1))  + 
  theme_bw() + ggtitle("WT Gran processed together") + 
  ylab("cell type") + 
  scale_fill_manual(values = wt_rad21_hue[2]) + 
  scale_color_manual(values = wt_rad21_hue[2])

rad21_gran_viol <- granulocytes.df %>% 
  dplyr::filter(genotype == "RAD21") %>% 
  dplyr::mutate(cluster = factor(seurat_clusters, levels = rad21_gran_order)) %>% 
  ggplot2::ggplot(., aes(cluster, monocle_pseudotime, fill = genotype)) + 
  geom_jitter(aes(cluster, monocle_pseudotime, col = genotype),  alpha=0.5, size =0.5, 
              position = position_jitterdodge()) + 
  geom_violin(scale="width", alpha = 0.5) + 
  geom_pointrange(mapping=aes(cluster, monocle_pseudotime), 
                  stat = "summary", 
                  fun.min = function(x) {quantile(x,0.25)}, 
                  fun.max = function(x) {quantile(x,0.75)}, 
                  fun="median", 
                  show.legend = F, size = 0.3, col = "darkblue", 
                  position=position_dodge(width=1))  + 
  theme_bw() + ggtitle("RAD21 Gran processed together") + 
  ylab("cell type") + 
  scale_fill_manual(values = wt_rad21_hue[1]) + 
  scale_color_manual(values = wt_rad21_hue[1])

pdf(paste0(working_dir, "gran_rad21_wt_pseudo_violin_2008.pdf"))
grid.arrange(wt_gran_viol, rad21_gran_viol)
dev.off()


### end of figure 4B

### Figure 4C: This figure was done using PRISM
### end of 4C

### Figure 4D: Granulocytes gene expression violin plots.

genes <- c( "lyz", "lect2l", "mmp13a.1", "adam8a",
           "spi1b", "cxcr4b", "cebpb", "cebpa")
clusters <- c(15,3,6,4)

gran.obj <- subset(geno.combined, 
                   subset = seurat_clusters %in% clusters & genotype %in% c("WT", "RAD21") )



gran_norm_exp <- GetAssayData(object = gran.obj,
                              assay = "RNA",
                              slot = "data")
# keeping only macrophages genes:

for (gene in genes ) {
  message(paste0("processing: ", gene))
  gene_norm_exp <- gran_norm_exp[gene,] %>% 
    as.data.frame() %>% 
    dplyr::rename(value = ".") %>% 
    dplyr::mutate(gene = gene) %>% 
    rownames_to_column("BC")
  
  norm_exp.df <- gran.obj@meta.data %>% 
    as.data.frame() %>%  
    rownames_to_column() %>% 
    dplyr::select(c(rowname, genotype, seurat_clusters)) %>%  
    dplyr::left_join(., gene_norm_exp , by = c("rowname" = "BC"))
  
  p <- norm_exp.df %>% mutate(cluster = factor(seurat_clusters, levels =  clusters)) %>% 
    ggplot(., aes(cluster, value, fill= genotype, col = genotype))  + 
    theme_bw() + 
    geom_violin(size = 0.5, scale = "width")  + 
    facet_grid(gene~.) + 
    geom_pointrange(mapping = aes(x = cluster, y = value), 
                    stat = "summary", 
                    fun.min = function(x) {quantile(x,0.25)}, 
                    fun.max = function(x) {quantile(x,0.75)}, 
                    fun="median", show.legend = F, size = 0.1, 
                    col = "darkblue", 
                    position = position_dodge(width = 0.9)) + 
    theme(strip.text.y.right = element_text(angle = 0))
  
  plot.pdf <- paste0(outdir,"gran_", gene, "_violinplot_1207.pdf")
  pdf(plot.pdf)
  print(p)
  dev.off()
}


#### Supplementary figures

### Supplementary figure 1A: Metascape pathway enrichment on Bulk RNAseq. 
### DEG with FDR <5%, split by regulation: up or down.
### end of Supplementary figure 1A

### Supplementary figure 1B: Expression level of key selected genes in WT non kidney clusters.

wt.nonkidney.obj <- subset(geno.combined, genotype == "WT" & iskidney == "N")
allcols <- c('0' = "#E16272", '1' = "#CF2B57", '2' = "#6DC7B4", 
             '3' = "#4D739F", '4'="#225FAC", '5'="#46B449", '6' = "#8FB2DA",
             '7' = "#F7941D", '8' = "#CD81B7", '9' = "#92278F", '10' = "#AC98C9", 
             '11' = "#9B9B9B", '12' = "#F7941D", '13' = "#9B9B9B", '14' = "#84C441", 
             '15' = "#3986C7", '16' =  "#FFF200", '17' = "#F7941D", '18' = "#9B9B9B", 
             '19' = "#9B9B9B", '20' = "#F7941D", '21'="#9B9B9B")

interesting_genes <- c("mki67", "pcna", "myca", "mpeg1.1", "grn2", "cpa5", "klf1", "znfl2a", "gata1a")

clusters_order <- c(16, 10, 8,9,0,1,15,3,6,4,7,12,17,20,2,5,14)
wt.nonkidney.norm.exp <- wt.nonkidney.obj[["RNA"]]@data %>% as.data.frame()

wt.nonkidney.norm.exp.interesting.df <- wt.nonkidney.norm.exp %>% rownames_to_column("genename") %>% dplyr::filter(genename %in% interesting_genes) %>% pivot_longer(!genename, names_to =  "BC") %>% dplyr::left_join(., wt.nonkidney.obj@meta.data %>%  rownames_to_column() %>% dplyr::select(c(rowname, seurat_clusters)) , by = c("BC" = "rowname")) %>% mutate(cluster = factor(seurat_clusters, levels = clusters_order))

pdf.file <- paste0(seuratdir, "figS1B_nonkidney_genes_panel_violin_WT.pdf")
pdf(pdf.file)
wt.nonkidney.norm.exp.interesting.df %>% 
  mutate(genes = factor(genename, levels = interesting_genes)) %>% 
  ggplot(., aes(cluster, value, fill = cluster)) + 
  geom_violin(scale = "width", trim = 0.2 )  + 
  theme_bw() + 
  scale_fill_manual(values = allcols) + 
  geom_pointrange(mapping = aes(x = cluster, y = value), 
                  stat = "summary", 
                  fun.min = function(x) {quantile(x,0.25)}, 
                  fun.max = function(x) {quantile(x,0.75)}, 
                  fun="median", show.legend = F, 
                  position = position_dodge(.200), size = 0.1) + 
  facet_grid(genes ~ . )

dev.off()




### end of Supplementary figure 1B

### Supplementary figure 1C: Done using PRISM: The percentage of the different maturing erythroblasts identified in the stained cytospins
### end of Supplementary figure 1C

### Supplementary figure 2: Heatmap for the top30 conserved markers by cluster.

# to draw heatmap we need to reorder the column according to what we want to plot
# we finally decided to produce only 1 heatmap with both nonkidney and kidney cells.
# make sure we are on RNA assay:

DefaultAssay(geno.combined)  <- "RNA"
# to build the heatmap we will use the conservedMarkers genes. Genes that are differentially expressed from 1 cluster compared to the rest.
# during the Seurat processing, we have exported a tab file per cluster.
# now we are reading it:

### Getting the conserved markers:
# need to loop on cluster and retrieve conservedMarkers data
# retrieving the cluster information:

conservedMarks.df <- data.frame()
for (cl in levels(Idents(geno.combined)) ) {
  message(paste0("processing cluster: ", cl, "<br>" ))
  comMarks.file <-paste0(seuratdir, "conservedMarkers0.5/conservedMarkers_", cl, "_res0.5.xls")
  message(paste0("Searching for: ", comMarks.file))
  if (file.exists(comMarks.file)) {
    comMarks.data <- read.table(comMarks.file,
                                sep = "\t",
                                stringsAsFactors = F,
                                header = T)
    comMarks.data$cluster <- cl
    conservedMarks.df <- rbind(conservedMarks.df, comMarks.data)
  }
  else {
    stop(paste0(" Missing: ", comMarks.file))
  }
}

# only getting top30 genes expressed by cluster based on the WT avg expression.
conserved_top30genes <- conservedMarks.df  %>% 
  group_by(cluster) %>% 
  top_n(n=30, wt = WT_avg_log2FC) %>% 
  distinct(rowname) %>% 
  pull(rowname)

# to these conserved top 30 genes we would like to show expression of known markers from litterature:
# we are loading the genes that we would like to see the gene names in the heatmap.

genes2show.df <- read.table(paste0(seuratdir, "gene2show.tsv"), 
                            sep = "\t", header = T, 
                            stringsAsFactors = F)

genes2show <- genes2show.df %>% 
  separate_rows(genes, sep = ", ") %>% 
  distinct(genes) %>% 
  pull(genes)

kidney2show.df <- read.table(paste0(seuratdir, "kidney2show.tsv"), 
                             sep = "\t", 
                             header = T, 
                             stringsAsFactors = F)


kidney2show <- kidney2show.df %>% 
  separate_rows(genes, sep = ", ") %>% 
  distinct(genes) %>% 
  pull(genes)

allgenes2show <- unique(c(genes2show, kidney2show))

# to plot the heatmap we need to have the scale.data of the genes we are planning to show expression.

#allgenes2show %in% VariableFeatures(geno.combined)

total_genes_whole_heatmap <- unique(c(conserved_top30genes, allgenes2show))
color2add <- ifelse(total_genes_whole_heatmap %in% allgenes2show, "black", "transparent")

rad21_wt.combined <-subset(geno.combined, genotype %in% c("RAD21", "WT"))

pdf(paste0(seuratdir,"allClusterHeatmapTop30conservedMarkers.pdf"))

DoHeatmap(rad21_wt.combined, features = total_genes_whole_heatmap, raster = F) + 
  theme(axis.text.y = element_text(color = rev(x = color2add))) + 
  NoLegend()

dev.off()

## this needs to be reordered. nonkidney cells first and then kidney at the end.
# to do so, we use a vector of Idents() to reorder the whole plot.
clusters_order <- as.character( c(0:10, 12, 14, 15, 16, 17, 20, 11, 13, 18, 19, 21))

Idents(rad21_wt.combined) <- clusters_order

pdf(paste0(seuratdir,"allClusterHeatmapTop30conservedMarkersSepKidney.pdf"))
DoHeatmap(rad21_wt.combined, features = total_genes_whole_heatmap, raster = F) + theme(axis.text.y = element_text(color = rev(x = color2add))) + NoLegend()
dev.off()

### end of supplementary figure 2

### Supplementary figure 3A: BP Gene Ontology for Erythrocyte

clusters <- c(  8, 9, 0, 1)
cutoff <- 0.5
sign_dges <- RNA_dges.df %>% 
  dplyr::filter(comparison == "Rad21_v_WT", p_val_adj <0.05, cluster %in% clusters ) %>%  
  mutate(reg = ifelse(avg_log2FC >0 , "up", "down"))  %>% 
  mutate(clusterF = factor(cluster, levels = clusters))

gene.df <- bitr(sign_dges$rowname,
                fromType = "SYMBOL",
                toType = "ENTREZID",
                OrgDb = OrgDB)

sign_dges <- sign_dges %>% 
  dplyr::left_join(., gene.df, by=c("rowname" = "SYMBOL"))

sign_dges.ck <- compareCluster(ENTREZID~clusterF+reg,
                               data=sign_dges,
                               fun="enrichGO", OrgDb = OrgDB,
                               ont = "BP",
                               readable = T)

path.xls <- paste0(pathwaydir,
                   "eryt_go_enrich_pathways_ck.xls")
sign_dges.ck %>% as.data.frame() %>% write.table(., path.xls,
                                                 quote = F,
                                                 sep = "\t")

sign_dges.ck.simp <- clusterProfiler::simplify(sign_dges.ck,
                                               cutoff=cutoff,
                                               by="p.adjust",
                                               select_fun = min)

### manual reordering of the pathways enrich.

path_ordered.file <- paste0(pathwaydir,
                            "ordering_eryt_pathways_1.tsv")
path_ordered <- read.table(path_ordered.file, sep = "\t",
                           header = T, stringsAsFactors = F)

dotplot.pdf <- paste0(pathwaydir,
                      "eryt_go_enrich_pathways_ck_simp05_min_ManualSelectionUpDown_1.pdf")
pdf(dotplot.pdf, width = 10)

sign_dges.ck.simp %>% 
  as.data.frame() %>% 
  dplyr::filter(Description %in% path_ordered$pathways) %>% 
  mutate(DescriptionF = factor(Description, levels = rev(path_ordered$pathways))) %>% 
  mutate(regF = factor(reg, levels = c("up", "down")))  %>% 
  mutate(clusterF = factor(clusterF, levels = clusters)) %>% 
  ggplot(., aes(clusterF, DescriptionF, col = -log10(p.adjust), size = Count)) + 
  geom_point() + scale_color_continuous(type = "viridis") + 
  theme_bw()  + 
  facet_grid(~regF) + 
  theme(axis.text.y.left = element_text(size = 6)) + 
  ggtitle("Eryt manually ordered 1  up/down")

dev.off()

#### end of Supplementary figure 3A

#### Supplementary figure 3B

clusters <- c( 15,3,6,4) 
cutoff <- 0.5
gran_sign_dges <- RNA_dges.df %>% 
  dplyr::filter(comparison == "Rad21_v_WT", p_val_adj <0.05, cluster %in% clusters ) %>%  
  mutate(reg = ifelse(avg_log2FC >0 , "up", "down")) %>% 
  mutate(clusterF = factor(cluster, levels = clusters))

gran_gene.df <- bitr(gran_sign_dges$rowname,
                     fromType = "SYMBOL",
                     toType = "ENTREZID",
                     OrgDb = OrgDB)

gran_sign_dges <- gran_sign_dges %>% 
  dplyr::left_join(., gran_gene.df, by=c("rowname" = "SYMBOL"))

gran_sign_dges.ck <- compareCluster(ENTREZID~clusterF+reg,
                                    data=gran_sign_dges,
                                    fun="enrichGO", OrgDb = OrgDB,
                                    ont = "BP",
                                    readable = T)

gran.xls <- paste0(pathwaydir, "granulocyte_go_enrich_pathways_ck.xls")

gran_sign_dges.ck %>% 
  as.data.frame() %>% 
  write.table(., gran.xls, sep = "\t", quote = F)

gran_sign_dges.ck.simp <- clusterProfiler::simplify(gran_sign_dges.ck,
                                                    cutoff=cutoff,
                                                    by="p.adjust",
                                                    select_fun = min)

### manual reordering of the pathways enrich.

gran_path_ordered.file <- paste0(pathwaydir,
                                 "ordering_gran_pathways_1.tsv")
gran_path_ordered <- read.table(gran_path_ordered.file, sep = "\t",
                                header = T, stringsAsFactors = F)

dotplot.pdf <- paste0(pathwaydir,
                      "gran_go_enrich_pathways_ck_simp05_min.pdf")
pdf(dotplot.pdf, width = 10)

gran_sign_dges.ck.simp %>% 
  as.data.frame() %>% 
  dplyr::filter(Description %in% gran_path_ordered$pathways) %>% 
  mutate(DescriptionF = factor(Description, levels = rev(gran_path_ordered$pathways))) %>% 
  mutate(regF = factor(reg, levels = c("up", "down")))  %>% 
  mutate(clusterF = factor(clusterF, levels = clusters)) %>% 
  ggplot(., aes(clusterF, DescriptionF, col = -log10(p.adjust), size = Count)) + 
  geom_point() + 
  scale_color_continuous(type = "viridis") + 
  theme_bw()  + 
  facet_grid(~regF) + 
  theme(axis.text.y.left = element_text(size = 6)) + 
  ggtitle("Granulocyte up/down")
dev.off()
