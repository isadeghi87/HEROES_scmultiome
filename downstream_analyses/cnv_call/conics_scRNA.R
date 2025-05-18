
setwd('~/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes.noindex/home/projects/pools/AG_Thongjuea')

library(biomaRt)
library(beanplot)
library(mixtools)
library(pheatmap)
library(zoo)
library(squash)
library(CONICSmat)
library(Rtsne)
library(scran)
library(infercnv)
library(Seurat)

# load data 
seurat_obj = readRDS("./Result/10x_multiomics/integrative_analysis/processed_obj/I070_032_seurat_annotated.rds") 

# Extract the raw count matrix from the Seurat object
raw_counts <- as.matrix(GetAssayData(seurat_obj, assay = "RNA", slot = "counts"))

# Get the cell metadata (e.g., cell type information)
cell_annotations <- data.frame(Cell = names(seurat_obj@active.ident),CellType = seurat_obj@active.ident)
colnames(cell_annotations) = c('Cell','CellType')
rownames(cell_annotations) = NULL


## gene ordering file
gene_location_data = read.delim('~/Downloads/hg38_gencode_v27.txt',header = F)

genes = intersect(rownames(raw_counts),gene_location_data[,1])

gene_loc = gene_location_data[match(genes,gene_location_data[,1]),]
colnames(gene_loc) = c("Gene",'Chr','Start','End')

## read chromosome arms  positions
regions=read.table("~/Downloads/chromosome_arm_positions_grch38.txt",sep="\t",row.names = 1,header = T)
head(regions,n=5)

##  get gene pos
gene_pos=getGenePositions(rownames(raw_counts))

## normalize
normFactor=calcNormFactors(raw_counts)
