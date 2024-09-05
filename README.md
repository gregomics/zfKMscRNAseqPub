# zfKMscRNAseqPub
Repository for scripts used for the publication of Zebrafish Kidney Marrow scRNA-seq analysis.


### Goal

The goal of this repository is to contain the scripts generated for the analysis
of the zebrafish kidney marrow transcriptome at both bulk and single cells levels
.

### Experiment design

Single cell transcriptomics (scRNA-Seq) from Kidney marrow zebrafish.

2 conditions (Rad21a+/-) and Wild Type samples. Each condition has two biological replicates (in fact a pool of 2-4 Adult Zebrafish Whole Kidney Marrow).

  * WT-1
  * WT-2
  * RAD21-1
  * RAD21-2

### scRNA-Seq preprocessing

Sequencing files (fastq) from 3′ GEM library and Gel bead Kit v3.1 (10x Genomics) were processed using CellRangeR and aligned against the Zebrafish transcriptome: GRCz11 version 98 only considering protein coding genes as well as lncRNA.

See the slurm script: scripts/map_using_cellRanger.sl

CellRangeR processing will create a directory for each sample. The most insteresting part is stored into outs/filtered_feature_bc_matrix
Where are stored 3 gzipped files:

  * barcodes.tsv.gz (list of cell barcodes)
  * matrix.mtx.gz (raw expression)
  * features.tsv.gz (geneID and gene names)

These files can be found on GEO  under the identifier: GSE275537

### From preprocessed files to cluster and annotation

Analysis will continue using R and Seurat package.
The subsequent R script:

```
R/zfkm_seurat_analysis.R
```

Will do the following steps:

  * create Seurat object
  * Quality control and filtering
  * Normalisation and finding variable genes (vst method)
  * Integrating the different genotypes using the approach published in https://doi.org/10.1016/j.cell.2019.05.031
  * on Integrated assay run the scaling, PCA and UMAP and create clusters with a resolution of 0.5.
  * get the conserved markers between clusters
  * get differentially expressed genes between genotypes in each cluster using DESeq2
  * annotate cluster using SCATCH and published data.

The main outcome of this R script is a RDS file containing the seurat object carrying all the information. This will be further used to plot the figures of the paper.  
