library(devtools)
install.packages("RcppArmadillo")
install.packages("RcppProgress")
install_github("linxihui/NNLM")

install_github("MacoskoLab/liger")

devtools::install_github("cole-trapnell-lab/garnett")
devtools::install_github("czythu/scCancer")

# Check the dependencies
suppressMessages(library(Seurat))
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(garnett))
suppressMessages(library(xgboost))
suppressMessages(library(ggplot2))
suppressMessages(library(ggsci))
suppressMessages(library(cowplot))
suppressMessages(library(viridis))
suppressMessages(library(magrittr))
suppressMessages(library(dplyr))
suppressMessages(library(edgeR))
library(devtools)
# Load all files in the folder
library(scCancer)


# Input data
sample_id = 'I007_072'
object = readRDS(paste0("./integrative_analysis/processed_obj/",sample_id,"_seurat_obj.rds"))

DefaultAssay(object) = 'RNA'
object <- NormalizeData(object = object,
                        normalization.method = "LogNormalize",
                        scale.factor = 10000,
                        verbose = F)
object <- FindVariableFeatures(object, selection.method = "vst", nfeatures = min(5000, length(rownames(object))), verbose = FALSE)
object <- ScaleData(object, vars.to.regress = c("nCount_RNA", "percent.mt"), verbose = FALSE)

### predict cell type
t.results <- predCellType(object@assays$RNA$data)
object$scCancer_pred = t.results$type.pred
p_cl = DimPlot(object, reduction = "umap.rna", label = TRUE, repel = TRUE, group.by = 'seurat_clusters',label.size = 5)   
sc_plot = DimPlot(object, reduction = "umap.rna", label = TRUE, repel = TRUE, group.by = 'scCancer_pred',label.size = 5)        

pp = p_cl | sc_plot



# Load model and prepare data
model.path <- paste0(system.file("txt", package = "scCancer"), "/sc_xgboost.model")
genes.path <- paste0(system.file("txt", package = "scCancer"), "/genes-scRNA-tcga-sorted.txt")
model.ref <- xgb.load(model.path)

# features <- read.table(genes.path)$V1
features <- as.list(read.table(genes.path))[[1]]
testdata <- t(as.matrix(object@assays$RNA$scale.data))

temp <- matrix(data = 0, nrow = nrow(testdata), ncol = length(features),
               dimnames = list(rownames(testdata), features))
current.features <- colnames(testdata)

for(j in 1:length(features)){
  if(features[j] %in% current.features){
    temp[,j] <- testdata[, features[j]]
  }
}

# Prediction
testdata <- xgb.DMatrix(testdata)
predict.label <- predict(model.ref, testdata)
predict.label[which(predict.label > 0.5)] <- "malignant"
predict.label[which(predict.label <= 0.5)] <- "nonMalignant"
table(predict.label)

# Visualization
object$malignant.label <- predict.label
object <- RunPCA(object, npcs = 30, verbose = FALSE)
object <- RunUMAP(object, reduction = "pca", dims = 1:30, verbose = FALSE)
DimPlot(object, group.by = "malignant.label")
