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
p.obj = readRDS('./integrative_analysis/processed_obj/pool14_filtered.Seurat.rds')
p='pool14'

## load demultiplex results
maj = read.delim(file ='./Demultiplexing_results/demuxafy/summary/RNA_ATAC_pool14_majority.tsv',
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

sample_id <- "I007_072"
p.sample <- subset(p.obj, subset = assignment == sample_id)

# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(p.sample) <- "ATAC"
p.sample <- RunTFIDF(p.sample)
p.sample <- FindTopFeatures(p.sample, min.cutoff = 'q0')
p.sample <- RunSVD(p.sample)
p.sample <- RunUMAP(p.sample, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

DimPlot(p.sample, reduction = "umap.atac", repel = TRUE) + ggtitle("ATAC")

#p.sample <- SCTransform(p.sample, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
# Perform standard analysis of each modality independently RNA analysis
DefaultAssay(p.sample) <- "RNA"

p.sample <- NormalizeData(p.sample)
p.sample <- FindVariableFeatures(p.sample)
p.sample <- ScaleData(p.sample)
p.sample <- RunPCA(p.sample)
p.sample <- RunUMAP(p.sample,dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

DimPlot(p.sample, reduction = "umap.rna", repel = TRUE) + ggtitle("RNA")

#####################################
p.sample <- FindMultiModalNeighbors(p.sample, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
p.sample <- RunUMAP(p.sample, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
p.sample <- FindClusters(p.sample, graph.name = "wsnn", algorithm = 1, verbose = FALSE,resolution = 0.8)

#########Plots######################
pl1 <- DimPlot(p.sample, reduction = "umap.rna", group.by = "seurat_clusters",label = TRUE, label.size = 5, repel = TRUE) + ggtitle("RNA")
pl2 <- DimPlot(p.sample, reduction = "umap.atac",group.by = "seurat_clusters", label = TRUE, label.size =5, repel = TRUE) + ggtitle("ATAC")
pl3 <- DimPlot(p.sample, reduction = "wnn.umap", group.by = "seurat_clusters", label = TRUE, label.size = 5, repel = TRUE) + ggtitle("WNN")

pp = pl1 + pl2 +pl3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))



# Normalize and find variable features
p.sample <- FindMultiModalNeighbors(p.sample, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
p.sample <- RunUMAP(p.sample, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
p.sample<- FindClusters(p.sample, graph.name = "wsnn", algorithm = 1, verbose = FALSE,resolution = 0.8)


#### find cell markers ####
#RNA
# For scRNAseq 
DefaultAssay(p.sample) <- "RNA"
cluster_markers_rna <- FindAllMarkers(p.sample, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Optionally, select top markers for each cluster for the heatmap
top_markers_rna <- cluster_markers_rna %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

# Re-scale data including the genes of interest
p.sample <- ScaleData(p.sample, features = top_markers_rna$gene)

# Create a heatmap for the top RNA markers
DoHeatmap(p.sample, features = top_markers_rna$gene) + NoLegend()

# For scATACseq
DefaultAssay(p.sample) <- "ATAC"
cluster_markers_atac <- FindAllMarkers(p.sample, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top_markers_atac <- cluster_markers_atac %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(p.sample, features = top_markers_atac$gene) + NoLegend()

## add pool number and sample
p.sample$pool = p
p.sample$sample = sample_id
#saveRDS(object = p.sample,paste0('./integrative_analysis/processed_obj/',p,'_majority_filtered.rds'))


DimPlot(p.sample, reduction = "wnn.umap", group.by = "seurat_clusters", 
        label = TRUE, label.size = 5, repel = TRUE)

# # ## find doublets using scrubletR 
# library(DoubletFinder)
# 
# 
# ## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
# homotypic.prop <- modelHomotypic(p.sample@meta.data$seurat_clusters)           ## ex: annotations <- p.sample@meta.data$ClusteringResults
# nExp_poi <- round(0.075*nrow(p.sample@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
# nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
# 
# ## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
# p.sample <- doubletFinder(p.sample, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
# p.sample <- doubletFinder(p.sample, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_391", sct = FALSE)
# p.pass = subset(p.sample,DF.classifications_0.25_0.09_193 =='Singlet')

## save the barcodes from this sample
saveRDS(p.sample,paste0("./integrative_analysis/processed_obj/",sample_id,"_seurat_obj.rds"))
barcodes= rownames(p.sample@meta.data)
}
