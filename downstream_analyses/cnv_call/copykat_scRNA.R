
library(infercnv)
library(copykat)

## load scmultiome data
setwd("~/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes.noindex/home/projects/pools/AG_Thongjuea/Result/10x_multiomics/cnv/copykat/I070_032/")

sample_id = 'I070_032'
seurat_obj = readRDS(paste0("./integrative_analysis/processed_obj/",sample_id,"_seurat_obj.rds"))

# Extract the raw count matrix from the Seurat object
raw_counts <- seurat_obj@assays$RNA$counts

copykat.test <- copykat(rawmat=raw_counts, 
                        id.type="S",
                        ngene.chr=5, 
                        win.size=25, KS.cut=0.1, 
                        sam.name=sample_id, 
                        distance="euclidean",
                        norm.cell.names="",
                        output.seg="TRUE", 
                        plot.genes="TRUE", 
                        genome="hg20",n.cores=4)

pred.test <- data.frame(copykat.test$prediction)
pred.test <- pred.test[which(pred.test$copykat.pred %in% c("aneuploid","diploid")),]  ##keep defined cells
CNA.test <- data.frame(copykat.test$CNAmat)

## heatmap 
my_palette <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 3, name = "RdBu")))(n = 999)

chr <- as.numeric(CNA.test$chrom) %% 2+1
rbPal1 <- colorRampPalette(c('black','grey'))
CHR <- rbPal1(2)[as.numeric(chr)]
chr1 <- cbind(CHR,CHR)

rbPal5 <- colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = "Dark2")[2:1])
com.preN <- pred.test$copykat.pred
pred <- rbPal5(2)[as.numeric(factor(com.preN))]

cells <- rbind(pred,pred)
col_breaks = c(seq(-1,-0.4,length=50),seq(-0.4,-0.2,length=150),seq(-0.2,0.2,length=600),seq(0.2,0.4,length=150),seq(0.4, 1,length=50))

heatmap.3(t(CNA.test[,4:ncol(CNA.test)]),dendrogram="r", distfun = function(x) parallelDist::parDist(x,threads =4, method = "euclidean"), hclustfun = function(x) hclust(x, method="ward.D2"),
          ColSideColors=chr1,RowSideColors=cells,Colv=NA, Rowv=TRUE,
          notecol="black",col=my_palette,breaks=col_breaks, key=TRUE,
          keysize=1, density.info="none", trace="none",
          cexRow=0.1,cexCol=0.1,cex.main=1,cex.lab=0.1,
          symm=F,symkey=F,symbreaks=T,cex=1, cex.main=4, margins=c(10,10))

legend("topright", paste("pred.",names(table(com.preN)),sep=""), pch=15,col=RColorBrewer::brewer.pal(n = 8, name = "Dark2")[2:1], cex=0.6, bty="n")


## check subpopulations of tumor cells
tumor.cells <- pred.test$cell.names[which(pred.test$copykat.pred=="aneuploid")]
tumor.mat <- CNA.test[, which(colnames(CNA.test) %in% tumor.cells)]
hcc <- hclust(parallelDist::parDist(t(tumor.mat),threads =4, method = "euclidean"), method = "ward.D2")
hc.umap <- cutree(hcc,2)

rbPal6 <- colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = "Dark2")[3:4])
subpop <- rbPal6(2)[as.numeric(factor(hc.umap))]
cells <- rbind(subpop,subpop)

heatmap.3(t(tumor.mat),dendrogram="r", distfun = function(x) parallelDist::parDist(x,threads =4, method = "euclidean"), hclustfun = function(x) hclust(x, method="ward.D2"),
          ColSideColors=chr1,RowSideColors=cells,Colv=NA, Rowv=TRUE,
          notecol="black",col=my_palette,breaks=col_breaks, key=TRUE,
          keysize=1, density.info="none", trace="none",
          cexRow=0.1,cexCol=0.1,cex.main=1,cex.lab=0.1,
          symm=F,symkey=F,symbreaks=T,cex=1, cex.main=4, margins=c(10,10))

legend("topright", c("c1","c2"), pch=15,col=RColorBrewer::brewer.pal(n = 8, name = "Dark2")[3:4], cex=0.9, bty='n')


## check with seurat metadata
cls = as.data.frame(seurat_obj$seurat_clusters)
cls$cell.names = rownames(cls)

merged.data = merge(cls,pred.test,by.x= 'cell.names')
