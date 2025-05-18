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
sample_id = 'I070_032'
sample_obj = readRDS("./integrative_analysis/processed_obj/I070_032_seurat_obj.rds")


### find cell markers ####
#RNA
# For scRNAseq 
DefaultAssay(sample_obj) <- "RNA"
cluster_markers_rna <- FindAllMarkers(sample_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Optionally, select top markers for each cluster for the heatmap
top_markers_rna <- cluster_markers_rna %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DimPlot(sample_obj, reduction = "wnn.umap",group.by = 'seurat_clusters',seed = 12, label = TRUE, label.size = 6) + NoLegend()+labs(title = '')


top_markers_rna$gene[top_markers_rna$cluster=='10']

for ( cl in unique(top_markers_rna$cluster)) {
  cat("Cluster:", cl, "\n")  # Print the cluster name
  print(top_markers_rna$gene[top_markers_rna$cluster == cl])  # Print genes for the current cluster
  cat("\n")  # Add a new line for better readability
}


# Integration of Markers, UMAP, and inferCNV:
#   Cluster 0:
#   Markers: "FAM189A1," "ABI3BP," "DKK1," "LINC00639," "ADGRV1," "C1orf61," "LINC01484," "CUX2"
# UMAP: Central and large cluster.
# CNV Pattern: Strong CNV signals.
# Interpretation: These markers are often involved in epithelial-mesenchymal transition (EMT), a process seen in tumor invasion. This cluster likely represents tumor epithelial cells or tumor-associated mesenchyme.
# Conclusion: Tumor-associated cluster.
# Cluster 1:
#   Markers: "NYAP2," "TGM3," "KLF5," "AQP3," "PLPPR5"
# UMAP: Fairly distinct from central tumor clusters.
# CNV Pattern: Minimal CNV alterations.
# Interpretation: KLF5 and AQP3 are associated with stem cell regulation and cellular differentiation, potentially pointing to normal immune or stromal cells.
# Conclusion: Non-tumor immune or stromal cells.
# Cluster 2:
#   Markers: "PCSK2," "GAD2," "GABRB1," "KCNN2," "NMNAT2"
# UMAP: Fairly separate.
# CNV Pattern: Minimal CNV signals.
# Interpretation: GAD2 and GABRB1 are linked with neural pathways, potentially indicative of neuronal-like stromal cells. No strong CNV signal suggests non-tumor status.
# Conclusion: Likely non-tumor stromal cells.
# Cluster 3:
#   Markers: "SYNPR," "CSMD3," "GRIN3A," "CTNNA2," "RAB27B"
# UMAP: Adjacent to central clusters.
# CNV Pattern: Moderate CNV.
# Interpretation: SYNPR and GRIN3A are neuronal markers, but CTNNA2 and RAB27B suggest involvement in cytoskeletal processes, which could also play a role in cancer cell invasion.
# Conclusion: Likely tumor-associated neural-like cells.
# Cluster 4:
#   Markers: "NEFL," "LVRN," "NEFM," "PAK5," "TACR3"
# UMAP: Clustered near tumor cells.
# CNV Pattern: Strong CNV signals.
# Interpretation: NEFL and NEFM are neuronal markers but often show abnormal expression in tumor contexts, particularly in sarcomas.
# Conclusion: Likely tumor-associated mesenchyme or neural-like tumor cells.
# Cluster 5:
#   Markers: "ROR2," "THBS2," "SFRP2," "POSTN," "TNFSF11"
# UMAP: Near central tumor clusters.
# CNV Pattern: Strong CNV signals.
# Interpretation: ROR2, THBS2, and POSTN are strongly associated with tumor mesenchymal/stroma and EMT processes.
# Conclusion: Likely tumor stroma or tumor mesenchyme.
# Cluster 6:
#   Markers: "HOXD10," "HOXD11," "SHISA6," "RYR1"
# UMAP: Somewhat isolated.
# CNV Pattern: Moderate CNV.
# Interpretation: HOXD10 and HOXD11 are associated with limb development but are also linked to tumor progression.
# Conclusion: Likely tumor mesenchyme.
# Cluster 7:
#   Markers: "CYBB," "MS4A4E," "ANKRD22," "RAB39A," "SIGLEC1"
# UMAP: Reference population.
# CNV Pattern: Minimal CNV.
# Interpretation: SIGLEC1 and MS4A4E are macrophage markers. This cluster represents macrophages used as a normal reference.
# Conclusion: Non-tumor (macrophages).
# Cluster 8:
#   Markers: "TSPEAR," "STAB2," "KRTAP10-4," "PROX1-AS1"
# UMAP: Moderately separated.
# CNV Pattern: Mild CNV changes.
# Interpretation: STAB2 and PROX1-AS1 point toward lymphatic endothelial cells.
# Conclusion: Likely non-tumor lymphatic endothelial cells.
# Cluster 9:
#   Markers: "EMCN," "VWF," "PLVAP," "RAMP3"
# UMAP: Distinct from other clusters.
# CNV Pattern: Minimal CNV.
# Interpretation: VWF and PLVAP are endothelial markers, indicating vascular endothelial cells.
# Conclusion: Likely non-tumor endothelial cells.
# Cluster 10:
#   Markers: "ESCO2," "KIF18B," "HJURP," "MYBL2," "KIF23"
# UMAP: Near tumor clusters.
# CNV Pattern: Strong CNV signals.
# Interpretation: These genes are linked with cell cycle progression and chromosomal instability, suggesting an aggressive tumor cell cluster.
# Conclusion: Tumor cell cluster.
# Cluster 11:
#   Markers: "THEMIS," "TNIP3," "CD8A," "SIRPG," "CD7"
# UMAP: Fairly distinct.
# CNV Pattern: Minimal CNV changes.
# Interpretation: CD8A and CD7 indicate cytotoxic T-cells.
# Conclusion: Non-tumor immune cells (T-cells).
# Cluster 12:
#   Markers: "TFAP2B," "SCHLAP1," "TDO2," "BBOX1-AS1"
# UMAP: Close to central tumor clusters.
# CNV Pattern: Strong CNV signals.
# Interpretation: TFAP2B and SCHLAP1 are associated with tumor progression.
# Conclusion: Tumor cell cluster.
# Cluster 13:
#   Markers: "LINC00989," "ABCC9," "STEAP4," "ADRA1B," "RGS5"
# UMAP: Moderately separated.
# CNV Pattern: Moderate CNV.
# Interpretation: RGS5 and ABCC9 are vascular markers, indicating pericytes or vascular smooth muscle cells.
# Conclusion: Likely non-tumor pericytes.
# Cluster 14:
#   Markers: "LINC02725," "GALR1," "SIGLEC15," "CALCR"
# UMAP: Isolated.
# CNV Pattern: Minimal CNV.
# Interpretation: SIGLEC15 and CALCR suggest smooth muscle cells or stroma.
# Conclusion: Likely non-tumor stroma or smooth muscle cells.
# Final Summary:
#   Tumor-Associated Clusters: 0, 3, 4, 5, 6, 10, 12
# Non-Tumor Clusters: 1, 2, 7, 8, 9, 11, 13, 14

## tumor normal
new_cluster_ids <- c(
  '0' = 'stroma',
  '1' = 'myoepithelial',
  '2' = 'neuroendocrine',
  '3' = 'neuronal',
  '4' = 'neuronal',
  '5' = 'tumor-associated stroma',
  '6' = 'tumor-associated mesenchyme',
  '7' = 'macrophage',
  '8' = 'endothelial',
  '9' = 'vascular endothelial',
  '10' = 'tumor',
  '11' = 'T-cell',
  '12' = 'Epithelial tumor',
  '13' = 'Vascular smooth muscle cells',
  '14' = 'smooth mucsle cell'
)

## annotate clusters
DefaultAssay(sample_obj) = 'RNA'
FeaturePlot(sample_obj, features = c("ESCO2",  "KIF18B", "DEPDC1",
                                     "KIF4A",  "HJURP",  "MYBL2",  "HMMR",   "ASF1B",  "KIF23",  "NEK2" ))
FeaturePlot(sample_obj, features = c("THBS2" ))



# Optionally, check if all clusters have been covered by the new labels
unique(Idents(sample_obj))

# Rename the clusters in the Seurat object using the new labels
sample_obj <- RenameIdents(sample_obj, new_cluster_ids)

# To visualize the changes, you can use DimPlot with the new identities
p_cl = DimPlot(sample_obj, reduction = "wnn.umap",group.by = 'seurat_clusters',seed = 12, label = TRUE, label.size = 6) + NoLegend()+labs(title = '')
p_cell = DimPlot(sample_obj, reduction = "wnn.umap", label = TRUE, label.size = 5) + NoLegend()
pp = p_cell / p_cl
ggsave(paste0("./integrative_analysis/figures/",sample_id,"_cell_annotation.pdf"),
       plot = pp,width = 8,height = 5.5)

## coverage plot

DefaultAssay(sample_obj) = 'ATAC'
pcov = CoveragePlot(sample_obj, region = "KIF4A",
             extend.downstream = 100, extend.upstream = 100,features = 'KIF4A',expression.assay = 'RNA')+
  theme_bw(base_family = 18)

ggsave("./EMseq/results/figures/report/coverage_plot.png",plot = pcov,width = 8,height = 5.5)

saveRDS(sample_obj,file = "./integrative_analysis/processed_obj/I070_032_seurat_annotated.rds")
