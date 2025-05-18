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
sample_obj = readRDS("./integrative_analysis/processed_obj/I007_072_seurat_obj.rds")

### find cell markers ####
#RNA
# For scRNAseq 
DefaultAssay(sample_obj) <- "RNA"
cluster_markers_rna <- FindAllMarkers(sample_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Optionally, select top markers for each cluster for the heatmap
top_markers_rna <- cluster_markers_rna %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)

sample_obj <- ScaleData(sample_obj, features = top_markers_rna$gene)
DimPlot(sample_obj, reduction = "umap.rna", group.by = "seurat_clusters",label = TRUE, label.size = 5, repel = TRUE) + ggtitle("RNA")



# Create a named vector with new cluster identities
new_cluster_ids <- c(
  '0' = 'tumor',
  '1' = 'myoepithelial',
  '2' = 'tumor',
  '3' = 'tumor',
  '4' = 'tumor',
  '5' = 'tumor',
  '6' = 'mesenchymes',
  '7' = 'macrophage',
  '8' = 'tumor',
  '9' = 'tumor',
  '10' = 'endothelial',
  '11' = 'lymphatic endothelial',
  '12' = 'unclear',
  '13' = 'epithelial',
  '17' = 'mesenchymes'
)

## annotate clusters
DefaultAssay(sample_obj) = 'RNA'
FeaturePlot(sample_obj, features = c("EWSR1",'MYRIP','EPHA6','ROR2','PDGFRA','MET','PAX1'))
FeaturePlot(sample_obj, features = c('HOXD10','HOXD11','PDGFRA','ATP10B','ATP10A','AC009652.1','ZNF385B','CDH18','SV2B','WASH3P'))
FeaturePlot(sample_obj, features = c('NKAIN3','RALYL','PRSS12','CAMK2D','SLIT2','ADGRL3','GRID2','TSNARE1','ERG','SMYD3'))
FeaturePlot(sample_obj, features = c('GRID2'))

# Optionally, check if all clusters have been covered by the new labels
unique(Idents(sample_obj))

# Rename the clusters in the Seurat object using the new labels
sample_obj <- RenameIdents(sample_obj, new_cluster_ids)

# To visualize the changes, you can use DimPlot with the new identities
DimPlot(sample_obj, reduction = "wnn.umap", label = TRUE, label.size = 6) + NoLegend()

cluster_markers_rna$gene[cluster_markers_rna$cluster=='1']
DefaultAssay(sample_obj) = 'ATAC'
pcov = CoveragePlot(sample_obj, region = "GRID2",idents = c("tumor", "endothelial", "macrophage"),
             extend.downstream = 100, extend.upstream = 100,features = 'GRID2',expression.assay = 'RNA')+
  theme_bw(base_family = 18)

ggsave("./EMseq/results/figures/report/coverage_plot.png",plot = pcov,width = 8,height = 5.5)

saveRDS(sample_obj,file = "./integrative_analysis/processed_obj/I070_032_seurat_annotated.rds")
