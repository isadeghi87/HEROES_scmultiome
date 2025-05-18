##################################
##################################
setwd('~/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes.noindex/home/projects/pools/AG_Thongjuea/Result/10x_multiomics/')

library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ggplot2)
library(dplyr)
library(stringr)
library(Matrix)
library(DoubletFinder)

##################################
#Load pools multi-omics RNA-Seq
dir = "~/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes.noindex/home/projects/pools/AG_Thongjuea/Result/10x_multiomics/CellRanger_ARC/"
pools = paste0('pool',14:22)

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
  
  p2[["percent.mt"]] <- PercentageFeatureSet(p2, pattern = "^MT-")# if this function gives error follow below code
  #mito_genes <- grep("^MT-", rownames(p2@assays$RNA@data), value = TRUE)
  #p2 <- PercentageFeatureSet(p2, features = mito_genes[mito_genes %in% rownames(p2)], col.name = "percent.mt",assay = 'RNA')
  
  # create ATAC assay and add it to the object
  p2[["ATAC"]] <- CreateChromatinAssay(
    counts = atac_counts,
    sep = c(":", "-"),
    fragments = frag.file,
    annotation = annotation
  )
  
  #################################
  # Increase the maximum allowed size for globals to 5 GiB
  options(future.globals.maxSize = 5 * 1024^3)
  
  DefaultAssay(p2) <- "ATAC"
  p2 <- NucleosomeSignal(p2)
  p2 <- TSSEnrichment(p2)
  
  vln.p = VlnPlot(
    object = p2,
    features = c("nCount_RNA","nCount_ATAC","percent.mt", "TSS.enrichment"),
    ncol = 4,group.by = "orig.ident",
    pt.size = 0.1,log = T)
  
  ggsave(filename = paste0('./integrative_analysis/figures/',p,'_vln_plot.pdf'), plot = vln.p,width = 10,height = 8.5)
  
  px1 <- FeatureScatter(p2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  px2 <- FeatureScatter(p2, feature1 = "nCount_ATAC", feature2 = "nFeature_ATAC")
  
  px1|px2
  
  # filter out low quality cells
  p2.pass <- subset(
    x = p2,
    subset = nCount_ATAC < 300000 &
      nCount_ATAC > 1000 &
      nCount_RNA < 300000 &
      nCount_RNA > 500 &
      TSS.enrichment > 0.8 &
      percent.mt < 15
  )

  nrow(p2.pass@meta.data)

  
# 
  # # ATAC analysis
  # # We exclude the first dimension as this is typically correlated with sequencing depth
  # DefaultAssay(p2.pass) <- "ATAC"
  # p2.pass <- RunTFIDF(p2.pass)
  # p2.pass <- FindTopFeatures(p2.pass, min.cutoff = 'q0')
  # p2.pass <- RunSVD(p2.pass)
  # p2.pass <- RunUMAP(p2.pass, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
  # 
  # DimPlot(p2.pass, reduction = "umap.atac", repel = TRUE) + ggtitle("ATAC")
  # 
  # #p2.pass <- SCTransform(p2.pass, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
  # # Perform standard analysis of each modality independently RNA analysis
  # DefaultAssay(p2.pass) <- "RNA"
  # 
  # p2.pass <- NormalizeData(p2.pass)
  # p2.pass <- FindVariableFeatures(p2.pass)
  # p2.pass <- ScaleData(p2.pass)
  # p2.pass <- RunPCA(p2.pass)
  # p2.pass <- RunUMAP(p2.pass,dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
  # 
  # DimPlot(p2.pass, reduction = "umap.rna", repel = TRUE) + ggtitle("RNA")
  # 
  # #####################################
  # p2.pass <- FindMultiModalNeighbors(p2.pass, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
  # p2.pass <- RunUMAP(p2.pass, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
  # p2.pass <- FindClusters(p2.pass, graph.name = "wsnn", algorithm = 1, verbose = FALSE,resolution = 0.8)
  # 
  # #########Plots######################
  # pl1 <- DimPlot(p2.pass, reduction = "umap.rna", group.by = "seurat_clusters",label = TRUE, label.size = 5, repel = TRUE) + ggtitle("RNA")
  # pl2 <- DimPlot(p2.pass, reduction = "umap.atac",group.by = "seurat_clusters", label = TRUE, label.size =5, repel = TRUE) + ggtitle("ATAC")
  # pl3 <- DimPlot(p2.pass, reduction = "wnn.umap", group.by = "seurat_clusters", label = TRUE, label.size = 5, repel = TRUE) + ggtitle("WNN")
  # 
  # pp = pl1 + pl2 +pl3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
  # 
  # ## find doublets using scrubletR 
  # # library(scrubletR)
  # # library(reticulate)
  # # ps.pass <- scrublet_R(seurat_obj = p2.pass)
  # ## pK Identification (ground-truth) ------------------------------------------------------------------------------------------
  # sweep.res <- paramSweep(p2.pass, PCs = 1:10, sct = FALSE)
  # gt.calls <- p2.pass@meta.data[rownames(sweep.res[[1]]), "GT"].   ## GT is a vector containing "Singlet" and "Doublet" calls recorded using sample multiplexing classification and/or in silico geneotyping results 
  # sweep.stats <- summarizeSweep(sweep.res, GT = TRUE, GT.calls = gt.calls)
  # bcmvn <- find.pK(sweep.stats)
  # 
  # 
  # ## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
  # homotypic.prop <- modelHomotypic(p2.pass@meta.data$seurat_clusters)           ## ex: annotations <- p2.pass@meta.data$ClusteringResults
  # nExp_poi <- round(0.075*nrow(p2.pass@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
  # nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  # 
  # ## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
  # p2.pass <- doubletFinder(p2.pass, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
  # p2.pass <- doubletFinder(p2.pass, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_391", sct = FALSE)
  # p2.pass = subset(p2.pass,DF.classifications_0.25_0.09_360=='Singlet')
  # 
  # 
  # ggsave(filename = paste0('./integrative_analysis/figures/',p,'_umap_clusters.pdf'),plot = pp,width = 10,height = 8.5)
  # 
  saveRDS(p2.pass,file =paste0('./integrative_analysis/processed_obj/',p,'_filtered.Seurat.rds'))
}
  