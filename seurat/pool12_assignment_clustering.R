setwd('~/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes.noindex/home/projects/pools/AG_Thongjuea/Result/10x_multiomics/')

library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ggplot2)
library(dplyr)
library(stringr)

## load seurat obj
if(TRUE){
p.obj = readRDS('./integrative_analysis/processed_obj/pool12_filtered.Seurat.rds')
p='pool12'

## load demultiplex results
maj = read.delim(file ='./Demultiplexing_results/demuxafy/summary/RNA_ATAC_pool12_majority.tsv',
                 header = T)
maj = maj[,c('barcode','majority_vote')]
colnames(maj)[1:2] = c("barcode",'assignment')

p.obj.meta = p.obj@meta.data
p.obj.meta$barcode = rownames(p.obj.meta)
p.obj.merged= merge(p.obj.meta,maj,all.x=T)

rownames(p.obj.merged) = p.obj.merged$barcode
p.obj.merged = p.obj.merged[,-c(1)]

p.obj@meta.data = p.obj.merged
#pl.maj <- DimPlot(p.obj, reduction = "wnn.umap", group.by = "assignment", label = T, label.size = 5, repel = TRUE) 

## remove doublet and unassigned
# Subset to remove doublets and unassigned cells
# p.obj <- subset(p.obj, subset = assignment != "doublet" & assignment != "unassigned")
sample_id <- "I070_032"
p.obj <- subset(p.obj, subset = assignment == sample_id)

# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(p.obj) <- "ATAC"
p.obj <- RunTFIDF(p.obj)
p.obj <- FindTopFeatures(p.obj, min.cutoff = 'q0')
p.obj <- RunSVD(p.obj)
p.obj <- RunUMAP(p.obj, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

DimPlot(p.obj, reduction = "umap.atac", repel = TRUE) + ggtitle("ATAC")

#p.obj <- SCTransform(p.obj, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
# Perform standard analysis of each modality independently RNA analysis
DefaultAssay(p.obj) <- "RNA"

p.obj <- NormalizeData(p.obj)
p.obj <- FindVariableFeatures(p.obj)
p.obj <- ScaleData(p.obj)
p.obj <- RunPCA(p.obj)
p.obj <- RunUMAP(p.obj,dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

DimPlot(p.obj, reduction = "umap.rna", repel = TRUE) + ggtitle("RNA")

#####################################
p.obj <- FindMultiModalNeighbors(p.obj, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
p.obj <- RunUMAP(p.obj, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
p.obj <- FindClusters(p.obj, graph.name = "wsnn", algorithm = 1, verbose = FALSE,resolution = 0.8)

#########Plots######################
pl1 <- DimPlot(p.obj, reduction = "umap.rna", group.by = "seurat_clusters",label = TRUE, label.size = 5, repel = TRUE) + ggtitle("RNA")
pl2 <- DimPlot(p.obj, reduction = "umap.atac",group.by = "seurat_clusters", label = TRUE, label.size =5, repel = TRUE) + ggtitle("ATAC")
pl3 <- DimPlot(p.obj, reduction = "wnn.umap", group.by = "seurat_clusters", label = TRUE, label.size = 5, repel = TRUE) + ggtitle("WNN")

pp = pl1 + pl2 +pl3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))



# Normalize and find variable features
p.obj <- FindMultiModalNeighbors(p.obj, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
p.obj <- RunUMAP(p.obj, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
p.obj<- FindClusters(p.obj, graph.name = "wsnn", algorithm = 1, verbose = FALSE,resolution = 0.8)


#### find cell markers ####
#RNA
# For scRNAseq 
DefaultAssay(p.obj) <- "RNA"
cluster_markers_rna <- FindAllMarkers(p.obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Optionally, select top m arkers for each cluster for the heatmap
top_markers_rna <- cluster_markers_rna %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

# Re-scale data including the genes of interest
p.obj <- ScaleData(p.obj, features = top_markers_rna$gene)

# Create a heatmap for the top RNA markers
DoHeatmap(p.obj, features = top_markers_rna$gene) + NoLegend()

# For scATACseq
DefaultAssay(p.obj) <- "ATAC"
cluster_markers_atac <- FindAllMarkers(p.obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top_markers_atac <- cluster_markers_atac %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(p.obj, features = top_markers_atac$gene) + NoLegend()

## add pool number and sample
p.obj$pool = p
p.obj$sample = sample_id
#saveRDS(object = p.obj,paste0('./integrative_analysis/processed_obj/',p,'_majority_filtered.rds'))


DimPlot(p.obj, reduction = "wnn.umap", group.by = "seurat_clusters", 
        label = TRUE, label.size = 5, repel = TRUE)

## save the barcodes from this sample
saveRDS(p.obj,"./integrative_analysis/processed_obj/I070_032_seurat_obj.rds")
barcodes= rownames(p.obj@meta.data)
write.table(barcodes,file = "./integrative_analysis/tables/I070_032_barcodes.tsv",quote = F,col.names = F,row.names = F)
}
