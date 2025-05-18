library(Seurat)
library(dplyr)
library(ggplot2)
library(infercnv)
library(GenomicRanges)

## load scmultiome data
setwd("~/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes.noindex/home/projects/pools/AG_Thongjuea/")
sample_id='I007_072'
seurat_obj = readRDS(paste0("./Result/10x_multiomics/integrative_analysis/processed_obj/",sample_id,"_seurat_annotated.rds")) 

# Extract the raw count matrix from the Seurat object
raw_counts <- GetAssayData(seurat_obj, assay = "ATAC", slot = "counts")

# Calculate gene activity matrix from scATAC-seq data using Signac
gene_activities <- GeneActivity(object = seurat_obj, assay = "ATAC")  # Assay is set to ATAC data

# Add the gene activity matrix to the Seurat object
seurat_obj[["gene_activity"]] <- CreateAssayObject(counts = gene_activities)

# Normalize the gene activity data
DefaultAssay(seurat_obj) <- "gene_activity"
seurat_obj <- NormalizeData(
  object = seurat_obj,
  assay = 'gene_activity',
  normalization.method = 'LogNormalize',
  scale.factor = median(seurat_obj@assays$gene_activity@counts)
)


#seurat_obj <- NormalizeData(seurat_obj)
#seurat_obj <- ScaleData(seurat_obj)


#METADATA VARIABLES
sample <- seurat_obj
DefaultAssay(sample) <- "gene_activity"

metacell_content <- 5 # Number of cells to be merged into metacells.


# METACELL GENERATION.

# Will store all the metacells. The test column will be removed at the end.
whole_metacells <- data.frame(test = rownames(sample), row.names = rownames(sample))


# Will store the complete annotation for the metacells.
whole_annotation <- data.frame(cluster_names = "test", row.names = "test")
meta_counter <- 0 # To keep a count of the metacells that are created.

# Generate a new metadata column storing the mapping cell-metacell.
sample[["metacell_mapping"]] <- "not_mapped"

# Here my sample has already been labelled, so the cluster labels are actually in the 'cluster_names' column. 
#Use 'seurat_clusters' instead if not labelled.

# Checking list:
#Idents(sample) <- sample$seurat_clusters

#sample@meta.data$seurat_clusters <- sample@meta.data
  

for (cluster_id in levels(sample@active.ident)){
  # custom version
  #for (cluster_id in levels(sample$)){
  print(sprintf("Computing metacells for cluster %s.", cluster_id))
  
  # Will store the metacells per cluster.
  metacells <- data.frame(test = rownames(sample), row.names = rownames(sample))
  #chunksample <- sample[, sample$cluster_names == cluster_id]
  # Subset the sample by each cluster ID.
  chunksample <- sample[, sample@active.ident == cluster_id]
  # Subset the sample by each cluster ID.
  #chunksample <- sample[, sample$seurat_clusters  == cluster_id]
  # Get the count data as a data frame and transpose it so columns are GENES and rows are CELLS.
  countdata <- t(as.data.frame(Seurat::GetAssayData(chunksample, slot = "counts")))
  
  # Get the possible amount of metacells.
  times <- trunc(dim(countdata)[1] / metacell_content)
  for (i in seq(1,times)){
    meta_counter <- meta_counter + 1
    # Generate slice points for each metacell. i.e: 1-5, 6-10, 11-15...
    start <- ((i -1) * metacell_content + 1)
    end <- i * metacell_content
    # Compute the slice as a data frame containing the sum of the subsetted cells. dims = 1 row (metacell), X columns (genes)
    slice <- as.data.frame(colSums(countdata[start:end, ]))
    # Get the name of the cells merged.
    cell_names <- rownames(countdata[start:end, ])
    # Add the metacell.
    col_name <- sprintf("metacell_%s", meta_counter)
    metacells[[col_name]] <- slice[,1]
    # Add the metacell mapping to the cells in the sample.
    sample$metacell_mapping[colnames(sample) %in% cell_names] <- col_name
  }
  # Delete the test column as we already have more than 1 column in our data frame.
  metacells[["test"]] <- NULL
  # Will contain the annotation of the generated metacells. Columns: cluster identities. Rows: each metacell.
  annotation <- data.frame(cluster_names = colnames(metacells), row.names = colnames(metacells))
  # Replace the dummy cluster_names column's values for the actual label for the cluster.
  annotation$cluster_names <- cluster_id
  # Add the annotation data and the metacell data to the "whole" dataframe. \
  # In the end: Number of Columns for metacell object = Number of rows for annotation object.
  whole_metacells <- cbind(whole_metacells, metacells)
  whole_annotation <- rbind(whole_annotation, annotation)
}

# Turn the names into characters for the sake of avoiding errors when subsetting.
# if (length(args) == 2) {
#   whole_annotation$cluster_names <- paste0("cl",as.character(whole_annotation$cluster_names))
# }

# Delete the test row from the global annotation data.
whole_annotation <- whole_annotation[!rownames(whole_annotation) %in% c("test"), , drop = FALSE]

# Delete the test column from the global metacell data.
whole_metacells$test <- NULL

# Path and name of the annotation file that will be used in the inferCNV call.
#cnv_analysis_folder <- "" # Path to the folder that will store the annotation file.
#dir.create(cnv_analysis_folder, recursive = TRUE)

annotation_file <- paste0("./Result/10x_multiomics/cnv/annotaions/",sample_id,"_cell_annotation_atac.txt")

# Save the annotation object.
utils::write.table(whole_annotation,
                   file = annotation_file,
                   sep = "\t",
                   row.names = TRUE,
                   col.names = FALSE,
                   quote = FALSE)

# Return the metacell object as a matrix (required for running inferCNV).
my_umi <- as.matrix(whole_metacells)

my.ref<-c('Endothelial','Macrophages')

infercnv_obj = CreateInfercnvObject(raw_counts_matrix=my_umi,
                                    annotations_file=annotation_file,
                                    delim="\t",
                                    gene_order_file='~/Downloads/hg38_gencode_v27.txt',
                                    ref_group_names=my.ref)

outdir = paste0("./Result/10x_multiomics/cnv/infercnv_atac_",sample_id)
# perform infercnv operations to reveal cnv signal
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir=outdir,  # dir is auto-created for storing outputs
                             cluster_by_groups=T,
                             smooth_method = "runmeans",
                             output_format = "pdf",
                             denoise=T,
                             HMM=FALSE, 
                             resume_mode = FALSE,
                             scale_data = TRUE
                             
)


