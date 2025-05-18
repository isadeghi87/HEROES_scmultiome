setwd('~/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes.noindex/home/projects/heroes-aya-pools/AG_Thongjuea/Result/10x_multiomics/')

library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ggplot2)
library(dplyr)
library(stringr)

## load seurat obj
p.obj = readRDS('./Seurat_objects/Pool3_filtered.Seurat.rds')

## load demultiplex results
demult = read.delim(file ='../../../temp_analysis/scRNA-Seq/Demultiplexing_results/scAvengers/results/clusters_doublet.tsv',
                 header = T)
# demult = subset(demult, status =='singlet')
#demult = demult[,c('barcode','assignment')]
demult$assignment[demult$status=='doublet']='doublet'
demult$assignment[demult$status=='unassigned']='unassigned'

meta = p.obj@meta.data
meta$barcode = rownames(meta)
merged= merge(meta,demult,all.x=T)

rownames(merged) = merged$barcode
merged = merged[,-c(1)]

p.obj@meta.data = merged

pl.assing <- DimPlot(p.obj, reduction = "wnn.umap", group.by = "assignment", label = T, label.size = 5, repel = TRUE) + 
  ggtitle("ATACseq_demultiplexing_pool3")

pl.status <- DimPlot(p.obj, reduction = "wnn.umap", group.by = "status", label = T, label.size = 5, repel = TRUE) + 
  ggtitle("ATACseq_demultiplexing_pool3")
mixPlot = pl.assing + pl.status

ggsave("./integrative_analysis/figures/pool3_atac_demultiplex.pdf",plot = mixPlot,
       width = 10,height = 8)

## remove doublet and unassigned
# Subset to remove doublets and unassigned cells
p.singlet <- subset(p.obj, subset = assignment != "doublet" & assignment != "unassigned")

# Normalize and find variable features
p.singlet <- FindMultiModalNeighbors(p.singlet, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
p.singlet <- RunUMAP(p.singlet, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
p.singlet<- FindClusters(p.singlet, graph.name = "wsnn", algorithm = 1, verbose = FALSE,resolution = 0.8)

p.singlet <- DimPlot(p.singlet, reduction = "wnn.umap", group.by = "assignment", 
                       label = TRUE, label.size = 5, repel = TRUE) + ggtitle("ATACseq_filtered_pool3")


ggsave("./integrative_analysis/figures/pool3_atac_demultiplex_filtered.pdf",plot = p.singlet,
       width = 10,height = 8)

DefaultAssay(p.singlet) <- "ATAC"
p.singlet <- RunTFIDF(p.singlet)
p.singlet <- FindTopFeatures(p.singlet, min.cutoff = 'q0')
p.singlet <- RunSVD(p.singlet)
p.singlet <- RunUMAP(p.singlet, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

atac = DimPlot(p.singlet, reduction = "umap.atac", repel = TRUE) + ggtitle("ATAC")


#p.singlet <- SCTransform(p.singlet, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
# Perform standard analysis of each modality independently RNA analysis
DefaultAssay(p.singlet) <- "RNA"

p.singlet <- NormalizeData(p.singlet)
p.singlet <- FindVariableFeatures(p.singlet)
p.singlet <- ScaleData(p.singlet)
p.singlet <- RunPCA(p.singlet)
p.singlet <- RunUMAP(p.singlet,dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

rna = DimPlot(p.singlet, reduction = "umap.rna", repel = TRUE) + ggtitle("RNA")

mixFilterd = atac+rna
ggsave("./integrative_analysis/figures/pool3_atac_demultiplex_filtered_mix.pdf",plot = mixFilterd,
       width = 10,height = 8)
