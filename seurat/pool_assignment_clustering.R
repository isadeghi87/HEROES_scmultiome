setwd('~/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes.noindex/home/projects/heroes-aya-pools/AG_Thongjuea/Result/10x_multiomics/')

library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ggplot2)
library(dplyr)
library(stringr)

## load seurat obj
p.obj = readRDS('./Seurat_objects/Pool4_filtered.Seurat.rds')
p.majority = p.vireo = p.souporcell = p.obj

## load demultiplex results
maj = read.delim(file = './Demultiplexing_results/downstream_analysis/multi_omix/results/pool4_1016717_majority.tsv',
                 header = T)
colnames(maj) = c("barcode",'assignment')

p.majority.meta = p.majority@meta.data
p.majority.meta$barcode = rownames(p.majority.meta)
p.majority.merged= merge(p.majority.meta,maj,all.x=T)

rownames(p.majority.merged) = p.majority.merged$barcode
p.majority.merged = p.majority.merged[,-c(1)]

p.majority@meta.data = p.majority.merged
pl.maj <- DimPlot(p.majority, reduction = "wnn.umap", group.by = "assignment", label = T, label.size = 5, repel = TRUE) + 
  ggtitle("majority votes- pool4")

## remove doublet and unassigned
# Subset to remove doublets and unassigned cells
p.majority <- subset(p.majority, subset = assignment != "doublet" & assignment != "unassigned")

# Normalize and find variable features
p.majority <- FindMultiModalNeighbors(p.majority, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
p.majority <- RunUMAP(p.majority, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
p.majority<- FindClusters(p.majority, graph.name = "wsnn", algorithm = 1, verbose = FALSE,resolution = 0.8)

pl.maj.filt <- DimPlot(p.majority, reduction = "wnn.umap", group.by = "assignment", 
               label = TRUE, label.size = 5, repel = TRUE) + ggtitle("majority-filtered-pool4")


###
vireo_rna.file = "./Demultiplexing_results/vireo/results/multiomix/RNA/WGS_BJ/q2000_dp100/pool4_1016717/donor_ids.tsv"
vireo_rna = read.table(file= vireo_rna.file,header=T)
vireo_rna<-vireo_rna[,c("cell","donor_id")]
colnames(vireo_rna)<-c("barcode","vireo_RNA_sample")

## rename sample IDs back to originals
#pool3
# vireo_rna$vireo_RNA_sample = str_replace_all(vireo_rna$vireo_RNA_sample,
#                                              c("AS-1058560-LR-69403" = "I010_013",
#                                                "AS-1058563-LR-69404" = "I005_014",
#                                                "AS-1058556-LR-69402" = "I023_037",
#                                                "AS-1058553-LR-69401" = "I020_001"))

##pool4
vireo_rna$vireo_RNA_sample = str_replace_all(vireo_rna$vireo_RNA_sample,
                                             c('AS-1058561-LR-69403' = 'I010_022', 
                                               'AS-1058553-LR-69401' = 'I020_004', 
                                               'AS-1058564-LR-69404' = 'I005_016',
                                               'AS-1058556-LR-69402' = 'I023_047'))




p.vireo.meta = p.vireo@meta.data
p.vireo.meta$barcode = rownames(p.vireo.meta)
p.vireo.meta.merged = merge(p.vireo.meta,vireo_rna,all.x=T)

rownames(p.vireo.meta.merged) = p.vireo.meta.merged$barcode
p.vireo.meta.merged = p.vireo.meta.merged[,-c(1)]
p.vireo@meta.data = p.vireo.meta.merged


pl.vireo <- DimPlot(p.vireo, reduction = "wnn.umap", group.by = "vireo_RNA_sample",
                    label = TRUE, label.size = 5, repel = TRUE) + ggtitle("vireo-pool4")

## remove doublet and unassigned
# Subset to remove doublets and unassigned cells
p.vireo <- subset(p.vireo, subset = vireo_RNA_sample != "doublet" & vireo_RNA_sample != "unassigned")

# Normalize and find variable features
p.vireo <- FindMultiModalNeighbors(p.vireo, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
p.vireo <- RunUMAP(p.vireo, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
p.vireo <- FindClusters(p.vireo, graph.name = "wsnn", algorithm = 1, verbose = FALSE,resolution = 0.8)

pl.vireo.filt <- DimPlot(p.vireo, reduction = "wnn.umap", group.by = "vireo_RNA_sample", 
               label = TRUE, label.size = 5, repel = TRUE) + ggtitle("vireo_filtered_pool4")



#### find cell markers ####
#RNA
# For scRNAseq 
DefaultAssay(p.majority) <- "RNA"
cluster_markers_rna <- FindAllMarkers(p.majority, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# Optionally, select top markers for each cluster for the heatmap
top_markers_rna <- cluster_markers_rna %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
# Re-scale data including the genes of interest
p.majority <- ScaleData(p.majority, features = top_markers_rna$gene)

# Create a heatmap for the top RNA markers
DoHeatmap(p.majority, features = top_markers_rna$gene) + NoLegend()

# For scATACseq
DefaultAssay(p.majority) <- "ATAC"
cluster_markers_atac <- FindAllMarkers(p.majority, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top_markers_atac <- cluster_markers_atac %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
DoHeatmap(p.majority, features = top_markers_atac$gene) + NoLegend()

## add pool number 
p.majority$pool = p
saveRDS(object = p.majority,paste0('./integrative_analysis/processed_obj/',p,'_majority_filtered.rds'))

## 
## to make the visualization easier, subset T cell clusters
cluster.names <- levels(p.majority)
cluster.names <- grep("0", cluster.names,value = TRUE)[1]
cluster <- subset(p.majority, idents = cluster.names)
CoveragePlot(cluster, region = 'PLCXD1', features = 'PLCXD1', assay = 'ATAC', expression.assay = 'RNA', peaks = FALSE)
