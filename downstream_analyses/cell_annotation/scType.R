## here we run scType to predict cell types
setwd('~/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes.noindex/home/projects/pools/AG_Thongjuea/Result/10x_multiomics/')


# load libraries and functions
lapply(c("dplyr","Seurat","HGNChelper","openxlsx"), library, character.only = T)
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R");
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

## load seurat obj

sample_obj = readRDS("./integrative_analysis/processed_obj/I070_032_seurat_obj.rds")

# DB file
db_ <- "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
tissue <- "Muscle" # e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus 

# prepare gene sets
gs_list <- gene_sets_prepare(db_, tissue)

##

# check Seurat object version (scRNA-seq matrix extracted differently in Seurat v4/v5)
seurat_package_v5 <- isFALSE('counts' %in% names(attributes(sample_obj[["RNA"]])));
print(sprintf("Seurat object %s is used", ifelse(seurat_package_v5, "v5", "v4")))

# extract scaled scRNA-seq matrix
DefaultAssay(sample_obj) = 'RNA'
scRNAseqData_scaled <- GetAssayData(sample_obj[['RNA']],slot = 'data')


# run ScType
es.max <- sctype_score(scRNAseqData = scRNAseqData_scaled, scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)

# NOTE: scRNAseqData parameter should correspond to your input scRNA-seq matrix. For raw (unscaled) count matrix set scaled = FALSE
# When using Seurat, we use "RNA" slot with 'scale.data' by default. Please change "RNA" to "SCT" for sctransform-normalized data,
# or to "integrated" for joint dataset analysis. To apply sctype with unscaled data, use e.g. sample[["RNA"]]$counts or sample[["RNA"]]@counts, with scaled set to FALSE.

# merge by cluster
cL_resutls <- do.call("rbind", lapply(unique(sample_obj@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(sample_obj@meta.data[sample_obj@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(sample_obj@meta.data$seurat_clusters==cl)), 10)
}))

sctype_scores <- cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] <- "Unknown"
print(sctype_scores[,1:3])


sample_obj@meta.data$sctype_classification = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  sample_obj@meta.data$sctype_classification[sample_obj@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

p_cl = DimPlot(sample_obj, reduction = "umap.rna", label = TRUE, repel = TRUE, group.by = 'seurat_clusters',label.size = 5)   
sc_plot = DimPlot(sample_obj, reduction = "umap.rna", label = TRUE, repel = TRUE, group.by = 'sctype_classification',label.size = 5)        

pp = p_cl | sc_plot


ggsave("./integrative_analysis/figures/I070_032_annotation_sctype.pdf",pp,width = 6,height = 4)
