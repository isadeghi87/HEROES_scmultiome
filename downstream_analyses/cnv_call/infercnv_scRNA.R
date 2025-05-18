
library(infercnv)

## load scmultiome data
setwd("~/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes.noindex/home/projects/pools/AG_Thongjuea/Result/10x_multiomics/")

sample_id = 'I007_072'
seurat_obj = readRDS(paste0("./integrative_analysis/processed_obj/",sample_id,"_seurat_obj.rds"))

# Extract the raw count matrix from the Seurat object
raw_counts <- GetAssayData(seurat_obj, assay = "RNA", slot = "counts")

# Get the cell metadata (e.g., cell type information)
cell_annotations <- data.frame(Cell = names(seurat_obj@active.ident),CellType = seurat_obj@active.ident)
colnames(cell_annotations) = c('Cell','CellType')
rownames(cell_annotations) = NULL
cell_annotations$CellType = paste0('cluster',cell_annotations$CellType)

# Save this annotation data to a file (this will be used by inferCNV)
annot_file=paste0("./cnv/annotaions/",sample_id,"_cell_annotation_file.txt")
write.table(cell_annotations,annot_file, sep="\t", quote=FALSE, row.names=FALSE, col.names=F)

## gene ordering file
gene_location_data = read.delim('~/Downloads/hg38_gencode_v27.txt',header = F)

genes = intersect(rownames(raw_counts),gene_location_data[,1])

gene_loc = gene_location_data[match(genes,gene_location_data[,1]),]
colnames(gene_loc) = c("Gene",'Chr','Start','End')

# Save the gene ordering file
write.table(gene_loc, "./cnv/annotaions/gene_ordering_file.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=F)

infercnv_obj <- CreateInfercnvObject(raw_counts_matrix = raw_counts,
                                     annotations_file = annot_file,
                                     delim = "\t",
                                     gene_order_file =  "./cnv/annotaions/gene_ordering_file.txt",
                                     ref_group_names = c('cluster5'))

infercnv_obj <- infercnv::run(infercnv_obj,
                              cutoff = 1,  # You can adjust this cutoff based on your dataset
                              out_dir = "cnv/infercnv_I007_072/",
                              cluster_by_groups = TRUE,
                              denoise = TRUE,
                              HMM = TRUE,
                              resume_mode = FALSE)  # Forces a fresh start


