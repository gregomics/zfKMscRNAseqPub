# module load R-bundle-Bioconductor
library(DropletUtils)
library(scDblFinder)
library(tidyverse)
resdir <- "/path_to_cellranger_directory/"
sampleIDs <- c("WT-1", "WT-2", "RAD21-1", "RAD21-2")

for (sampleID in sampleIDs) {
  dir <- paste0(resdir, sampleID, 
		"/outs/filtered_feature_bc_matrix/")
  if (dir.exists(dir)) {
    cat(paste0("found ", dir, "\n"))
    sce <- read10xCounts(dir)
    sce <- scDblFinder(sce)
    ### exporting table
    file <- paste0(resdir, sampleID, 
		   "_scDblFinder.tsv")
    colData(sce) %>% as.data.frame() %>% dplyr::select(Barcode, scDblFinder.class) %>% write.table(., file, quote = F)
    singlets <- colData(sce) %>% as.data.frame() %>% dplyr::filter(scDblFinder.class == "singlet") %>% nrow()
    doublets <- colData(sce) %>% as.data.frame() %>% dplyr::filter(scDblFinder.class == "doublet") %>% nrow()
    ## plotting the the number of barcide per scDblFinder class (e.g. singlet or doublet)
    rate <- (doublets /(doublets + singlets)) *100
    title <- paste0(sampleID, ": doublet rate: ", rate)
    q <- colData(sce) %>% as.data.frame() %>% dplyr::select(Barcode, scDblFinder.class) %>% group_by(scDblFinder.class) %>% summarise(count = n() ) %>% ggplot(., aes(scDblFinder.class, count, fill = scDblFinder.class)) + geom_bar(stat = "identity") + theme_bw() + ggtitle(title)
    pdf.file <- paste0(resdir, sampleID, 
		       "_scDblFinder.pdf")
    pdf(pdf.file)
    print(q)
    dev.off()
  } else {
     cat(paste0("dir: ", dir, " is missing!"))
  }
}
