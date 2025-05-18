##################################
##################################
setwd('~/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes.noindex/home/projects/heroes-aya-pools/AG_Thongjuea/Result/10x_multiomics/')

library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ggplot2)
library(dplyr)
library(stringr)
library(Matrix)

##################################
#Load pools multi-omics RNA-Seq
pools = paste0('pool',5:8)
dir = "~/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes.noindex/home/projects/heroes-aya-pools/AG_Thongjuea/Code/10x_multiomics/Sadeghi/cellRanger_arc/"

for (p in pools){
  h5 = list.files(paste0(dir,p),pattern = 'filtered_feature_bc_matrix.h5',recursive = T,full.names = T)
  h5.p = h5[grep(p,h5)]
  inputdata.10x <- Read10X_h5(h5)
  
  # frag.file <- paste0("./CellRanger_arc/1016713_",p,"_hg38/outs/atac_fragments.tsv.gz")
  frag = list.files(paste0(dir,p),pattern = 'atac_fragments.tsv.gz',recursive = T,full.names = T)
  frag.file = frag[1]
  
  # extract RNA and ATAC data
  rna_counts <- inputdata.10x$`Gene Expression`
  atac_counts <- inputdata.10x$Peaks
  
  ##################################
  # get gene annotations for hg38
  annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
  seqlevels(annotation) <- paste0('chr', seqlevels(annotation))
  
  # create a Seurat object containing the RNA data
  p2 <- CreateSeuratObject(
    counts = rna_counts,
    assay = "RNA"
  )
  
  #p2[["percent.mt"]] <- PercentageFeatureSet(p2, pattern = "^MT-")# if this function gives error follow below code
  mito_genes <- grep("^MT-", rownames(p2@assays$RNA@data), value = TRUE)
  p2 <- PercentageFeatureSet(p2, features = mito_genes[mito_genes %in% rownames(p2)], col.name = "percent.mt",assay = 'RNA')
  
  # create ATAC assay and add it to the object
  p2[["ATAC"]] <- CreateChromatinAssay(
    counts = atac_counts,
    sep = c(":", "-"),
    fragments = frag.file,
    annotation = annotation
  )
  
  #################################
  
  DefaultAssay(p2) <- "ATAC"
  p2 <- NucleosomeSignal(p2)
  p2 <- TSSEnrichment(p2)
  
  vln.p = VlnPlot(
    object = p2,
    features = c("nCount_RNA","nCount_ATAC","percent.mt", "TSS.enrichment"),
    ncol = 4,group.by = "orig.ident",
    pt.size = 0.1,log = T
  )
  ggsave(filename = paste0('./integrative_analysis/figures/',p,'_vln_plot.pdf'),plot = vln.p,width = 10,height = 8.5)
  
  px1 <- FeatureScatter(p2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  px2 <- FeatureScatter(p2, feature1 = "nCount_ATAC", feature2 = "nFeature_ATAC")
  
  
  # filter out low quality cells
  p2.pass <- subset(
    x = p2,
    subset = nCount_ATAC < 300000 &
      nCount_RNA < 150000 &
      nCount_ATAC > 5000 &
      nCount_RNA > 1000 &
      TSS.enrichment > 1 &
      percent.mt < 20
  )
  
  
  
  # call peaks using MACS2
  #peaks <- CallPeaks(p2)
  
  # remove peaks on nonstandard chromosomes and in genomic blacklist regions
  #peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
  #peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)
  
  # quantify counts in each peak
  #macs2_counts <- FeatureMatrix(
  #  fragments = Fragments(pbmc),
  #  features = peaks,
  #  cells = colnames(pbmc)
  #)
  
  # create a new assay using the MACS2 peak set and add it to the Seurat object
  #pbmc[["peaks"]] <- CreateChromatinAssay(
  #  counts = macs2_counts,
  #  fragments = fragpath,
  #  annotation = annotation
  #)
  
  # ATAC analysis
  # We exclude the first dimension as this is typically correlated with sequencing depth
  DefaultAssay(p2.pass) <- "ATAC"
  p2.pass <- RunTFIDF(p2.pass)
  p2.pass <- FindTopFeatures(p2.pass, min.cutoff = 'q0')
  p2.pass <- RunSVD(p2.pass)
  p2.pass <- RunUMAP(p2.pass, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
  
  DimPlot(p2.pass, reduction = "umap.atac", repel = TRUE) + ggtitle("ATAC")
  
  
  #p2.pass <- SCTransform(p2.pass, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
  # Perform standard analysis of each modality independently RNA analysis
  DefaultAssay(p2.pass) <- "RNA"
  
  p2.pass <- NormalizeData(p2.pass)
  p2.pass <- FindVariableFeatures(p2.pass)
  p2.pass <- ScaleData(p2.pass)
  p2.pass <- RunPCA(p2.pass)
  p2.pass <- RunUMAP(p2.pass,dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
  
  DimPlot(p2.pass, reduction = "umap.rna", repel = TRUE) + ggtitle("RNA")
  #####################################
  p2.pass <- FindMultiModalNeighbors(p2.pass, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
  p2.pass <- RunUMAP(p2.pass, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
  p2.pass <- FindClusters(p2.pass, graph.name = "wsnn", algorithm = 1, verbose = FALSE,resolution = 0.8)
  
  #########Plots######################
  pl1 <- DimPlot(p2.pass, reduction = "umap.rna", group.by = "seurat_clusters",label = TRUE, label.size = 5, repel = TRUE) + ggtitle("RNA")
  pl2 <- DimPlot(p2.pass, reduction = "umap.atac",group.by = "seurat_clusters", label = TRUE, label.size =5, repel = TRUE) + ggtitle("ATAC")
  pl3 <- DimPlot(p2.pass, reduction = "wnn.umap", group.by = "seurat_clusters", label = TRUE, label.size = 5, repel = TRUE) + ggtitle("WNN")
  
  pp = pl1 + pl2 +pl3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))

  ggsave(filename = paste0('./integrative_analysis/figures/',p,'_umap_clusters.pdf'),plot = pp,width = 10,height = 8.5)
  
  saveRDS(p2.pass,file =paste0('./integrative_analysis/processed_obj/',p,'_filtered.Seurat.rds'))
}
  
  ########Vireo RNA demultiplexing####
  vireo_rna.file = "./Demultiplexing_results/vireo/results/multiomix/RNA/WGS_BJ/q2000_dp100/pool3_1016715/donor_ids.tsv"
  vireo_rna = read.table(file= vireo_rna.file,header=T)
  vireo_rna<-vireo_rna[,c("cell","donor_id")]
  colnames(vireo_rna)<-c("barcode","vireo_RNA_sample")
  
  ## rename sample IDs back to originals
  vireo_rna$vireo_RNA_sample = str_replace_all(vireo_rna$vireo_RNA_sample,
                                               c("AS-1058560-LR-69403" = "I010_013",
                                                 "AS-1058563-LR-69404" = "I005_014",
                                                 "AS-1058556-LR-69402" = "I023_037",
                                                 "AS-1058553-LR-69401" = "I020_001"))
  
  # ########Vireo ATAC demultiplexing####
  # vireo_atac.file<-"./Demultiplexing_results/vireo/results/multiomix/RNA/WGS_BJ/q2000_dp100/pool3_1016715/donor_ids.tsv"
  # vireo_atac<-read.table(file=vireo_atac.file,header=T)
  # vireo_atac<-vireo_atac[,c("barcode","status","assignment")]
  # colnames(vireo_atac)<-c("barcode","vireo_ATAC_doublet_status","vireo_ATAC_sample")
  
  
  p2.meta = p2.pass@meta.data
  p2.meta$barcode = rownames(p2.meta)
  p2.meta.merged = merge(p2.meta,vireo_rna,all.x=T)
  # p2.meta.merged<-merge(p2.meta.merged,vireo_atac,all.x=T)
  
  rownames(p2.meta.merged) = p2.meta.merged$barcode
  p2.meta.merged = p2.meta.merged[,-c(1)]
  
  p2.f<-p2.pass
  p2.f@meta.data<-p2.meta.merged
  
  
  # pl1 <- DimPlot(p2.pass, reduction = "umap.rna", group.by = "seurat_clusters",label = TRUE, label.size = 5, repel = TRUE) + ggtitle("RNA")
  # pl2 <- DimPlot(p2.pass, reduction = "umap.atac",group.by = "seurat_clusters", label = TRUE, label.size =5, repel = TRUE) + ggtitle("ATAC")
  
  pl5 <- DimPlot(p2.f, reduction = "wnn.umap", group.by = "vireo_RNA_sample", label = FALSE, label.size = 5, repel = TRUE) + ggtitle("")
  # pl6 <- DimPlot(p2.f, reduction = "wnn.umap", group.by = "vireo_ATAC_sample", label = FALSE, label.size = 5, repel = TRUE) + ggtitle("")
  
  pl5 + pl6 & theme(plot.title = element_text(hjust = 0.5))
  ################################
  #p2.f@meta.data$vireo_RNA_sample_1 <- p2.f@meta.data$vireo_RNA_sample
  # p2.f@meta.data$vireo_RNA_sample_1[p2.f@meta.data$vireo_RNA_doublet_status=="doublet"]<-"doublet"
  # p2.f@meta.data$vireo_RNA_sample_1[p2.f@meta.data$vireo_RNA_doublet_status=="unassigned"]<-"unassigned"
  
  pl7 <- DimPlot(p2.f, reduction = "wnn.umap", group.by = "vireo_RNA_sample", 
                 label = FALSE, label.size = 5, repel = TRUE) + ggtitle("")
  
  
  saveRDS(p2.meta.merged,file="./../../Result/10x_multiomics/Seurat_objects/pool3_merged.rds")
  
  
  ## remove doublet and unassigned
  # Subset to remove doublets and unassigned cells
  p2.f <- subset(p2.f, subset = vireo_RNA_sample != "doublet" & vireo_RNA_sample != "unassigned")
  
  
  # Normalize and find variable features
  p2.f <- FindMultiModalNeighbors(p2.f, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
  p2.f <- RunUMAP(p2.f, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
  p2.f<- FindClusters(p2.f, graph.name = "wsnn", algorithm = 1, verbose = FALSE,resolution = 0.8)
  
  pl8 <- DimPlot(p2.f, reduction = "wnn.umap", group.by = "vireo_RNA_sample", 
                 label = FALSE, label.size = 5, repel = TRUE) + ggtitle("")
  #### find cell markers ####
  #RNA
  # For scRNAseq 
  DefaultAssay(p2.f) <- "RNA"
  cluster_markers_rna <- FindAllMarkers(p2.f, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  # Optionally, select top markers for each cluster for the heatmap
  top_markers_rna <- cluster_markers_rna %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
  # Re-scale data including the genes of interest
  p2.f <- ScaleData(p2.f, features = top_markers_rna$gene)
  
  # Create a heatmap for the top RNA markers
  DoHeatmap(p2.f, features = top_markers_rna$gene) + NoLegend()
  
  # For scATACseq
  DefaultAssay(p2.f) <- "ATAC"
  cluster_markers_atac <- FindAllMarkers(p2.f, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  top_markers_atac <- cluster_markers_atac %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
  DoHeatmap(p2.f, features = top_markers_atac$gene) + NoLegend()
  
  # cluster_markers_atac_clust <- FindAllMarkers(p2.f, ident.1 = cluster1, ident.2 = cluster2) # Specify clusters
  
  # Visualization of top markers
  VlnPlot(p2.f, features = head(rownames(cluster_markers_rna), 5))
  VlnPlot(p2.f, features = head(rownames(cluster_markers_atac), 5))
  
  
  ### Integration of scRNAseq and scATACseq Data ####
  # Assuming `p2.f` is your Seurat object with both RNA and ATAC assays
  # Perform SCTransform normalization
  p2.f <- SCTransform(p2.f, verbose = FALSE)
  
  # Select integration features
  p2.f <- SelectIntegrationFeatures(object = listp2.f, nfeatures = 3000)
  
  # Preparing for SCT integration
  p2.f <- PrepSCTIntegration(object = p2.f, anchor.features = rna.features)
  
  DefaultAssay(p2.f) <- "RNA"
  rna.features <- SelectIntegrationFeatures(object.list = list(p2.f), nfeatures = 3000)
  p2.integ <- PrepSCTIntegration(object = list(p2.f), anchor.features = rna.features)
  label = FALSE, label.size = 5, repel = TRUE) + ggtitle("")

# Run the PCA, UMAP and clustering on integrated data
p2.f <- ScaleData(p2.f)
p2.f <- RunPCA(p2.f)
p2.f <- RunUMAP(p2.f, reduction = "pca", dims = 1:30)
p2.f <- FindNeighbors(p2.f, dims = 1:30)
p2.f <- FindClusters(p2.f)

saveRDS(object = p2.f,file = './Demultiplexing_results/downstream_analysis/multi_omix/results/pool3_integrated.rds')

### trajectory analysis ####
library(monocle3)

# Convert Seurat object to Monocle's CellDataSet
# Convert Seurat object to SingleCellExperiment
p2.f_monocle <- as.SingleCellExperiment(p2.f)

# Preprocess the data
p2.f_monocle <- monocle3::preprocess_cds(p2.f_monocle, num_dim = 100)
p2.f_monocle <- reduce_dimension(p2.f_monocle)

# Order cells in pseudotime
cds <- orderCells(cds)

# Plotting the trajectory
plot_cell_trajectory(cds)

library(slingshot)

# Extracting required data from Seurat
sce <- as.SingleCellExperiment(p2.f)

# Running Slingshot for trajectory inference
sce <- slingshot(sce, clusterLabels = 'seurat_clusters', reducedDim = 'UMAP')

# Plotting the trajectory
plot(reducedDims(sce)$UMAP, col = rainbow(length(unique(sce$seurat_clusters)))[factor(sce$seurat_clusters)])
lines(SlingshotDataSet(sce), lwd = 2, type = 'lineages')

}


