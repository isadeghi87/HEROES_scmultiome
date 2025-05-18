## annotate cell types for sample "I070_032" from pool12
setwd('~/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes.noindex/home/projects/pools/AG_Thongjuea/Result/10x_multiomics/')

library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ggplot2)
library(dplyr)
library(stringr)

## read seurat obj
sample_id = 'I007_072'
sample_obj = readRDS(paste0("./integrative_analysis/processed_obj/",sample_id,"_seurat_obj.rds"))

### find cell markers ####
#RNA
# For scRNAseq 
DefaultAssay(sample_obj) <- "RNA"
cluster_markers_rna <- FindAllMarkers(sample_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Optionally, select top markers for each cluster for the heatmap
top_markers_rna <- cluster_markers_rna %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)

sample_obj <- ScaleData(sample_obj, features = top_markers_rna$gene)
DimPlot(sample_obj, reduction = "umap.rna", group.by = "seurat_clusters",label = TRUE, label.size = 5, repel = TRUE) + ggtitle("RNA")


for (cl in unique(top_markers_rna$cluster)) {
  cat("Cluster:", cl, "\n")  # Print the cluster name
  print(top_markers_rna$gene[top_markers_rna$cluster == cl])  # Print genes for the current cluster
  cat("\n")  # Add a new line for better readability
}


##gpt prediction
# Cluster 0: Mesenchymal or Fibroblast-like Cells
# Cluster 1: Epithelial or Neuroendocrine Cells
# Cluster 2: Proliferating Cells or Stromal Cells
# Cluster 3: Endothelial or Mesenchymal Cells
# Cluster 4: Proliferating Tumor Cells
# Cluster 5: Endothelial Cells or Myeloid Cells (Macrophages)

# Create a named vector with new cluster identities
new_cluster_ids <- c(
  '0' = 'Mesenchyme',
  '1' = 'Epithelial',
  '2' = 'Proliferating stromal Cells',
  '3' = 'Endothelial',
  '4' = 'Proliferating tumor Cells',
  '5' = 'Macrophages'
)

cluster_markers_rna$gene[cluster_markers_rna$cluster=='1']
## annotate clusters
DefaultAssay(sample_obj) = 'RNA'
FeaturePlot(sample_obj, features = c('KIF18B','MELK','FOXO1','PAX3','PAX7'))
FeaturePlot(sample_obj, features = c('PTEN'))

# Optionally, check if all clusters have been covered by the new labels
unique(Idents(sample_obj))

# Rename the clusters in the Seurat object using the new labels
sample_obj <- RenameIdents(sample_obj, new_cluster_ids)

# To visualize the changes, you can use DimPlot with the new identities
p_cl = DimPlot(sample_obj, reduction = "wnn.umap",group.by = 'seurat_clusters',seed = 12, label = TRUE, label.size = 6) + NoLegend()+labs(title = '')
p_cell = DimPlot(sample_obj, reduction = "wnn.umap", label = TRUE, label.size = 6) + NoLegend()
pp = p_cl | p_cell

ggsave(paste0("./integrative_analysis/figures/",sample_id,"_cell_annotation.pdf"),
              plot = pp,width = 8,height = 5.5)

DefaultAssay(sample_obj) = 'ATAC'
pcov = CoveragePlot(sample_obj, region = "MELK",
             extend.downstream = 100, extend.upstream = 100,features = 'MELK',expression.assay = 'RNA')+
  theme_bw(base_family = 18)
ggsave(paste0("./integrative_analysis/figures/",sample_id,"_MELK_coverage_plot.pdf"),
       plot = pcov,width = 8,height = 5.5)

(pcov2 = CoveragePlot(sample_obj, region = "PAX7",
                    extend.downstream = 100, extend.upstream = 100,features = 'PAX7',expression.assay = 'RNA')+
  theme_bw(base_family = 18))

ggsave(paste0("./integrative_analysis/figures/",sample_id,"_KIF18B_coverage_plot.pdf"),
       plot = pcov2,width = 8,height = 5.5)

saveRDS(sample_obj,paste0("./integrative_analysis/processed_obj/",sample_id,"_seurat_annotated.rds"))
