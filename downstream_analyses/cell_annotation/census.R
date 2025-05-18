setwd('~/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes.noindex/home/projects/pools/AG_Thongjuea/Result/10x_multiomics/')

library(Census)
library(dplyr)
library(Seurat)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)


## fix 
predict_node_M = function(obj, model, node, get_prob = T, allowed_nodes = NULL){
  
  idx = which(model$hierarchy_mat == node, arr.ind = T)
  new_ids = model$hierarchy_mat[idx[1,1] + 1, idx[,2]] %>% as.character() %>% unique()
  
  if(!is.null(allowed_nodes) & any(new_ids %in% allowed_nodes == F)){
    obj = obj[, Idents(obj) == node]
    pred_res = data.frame(barcode = colnames(obj),
                          pred = new_ids[new_ids %in% allowed_nodes],
                          prob = NA,
                          census_clusters = obj$census_clusters,
                          u1 = obj@reductions$umap@cell.embeddings[,1],
                          u2 = obj@reductions$umap@cell.embeddings[,2])
  }
  
  if(any(new_ids %in% allowed_nodes == F) & get_prob == F){
    return(pred_res)
  } else {
    g = intersect(model$markers[[node]]$gene, rownames(obj))
    temp_seurat = obj[g, Idents(obj) == node]
    
    x2 = temp_seurat@assays$RNA$counts %>% as.matrix()
    
    # add missing genes
    g = setdiff(model$markers[[node]]$gene, rownames(obj))
    if(length(g) > 0){
      x2 = rbind(x2, matrix(0, nrow = length(g), ncol = ncol(x2)))
      rownames(x2) = c(rownames(temp_seurat), g)
    }
    
    x2 = x2[model$markers[[node]]$gene, ]
    
    x2[x2 == 0] = NA
    
    x2 = apply(x2, 2, dplyr::ntile, 100) %>% t()
    colnames(x2) = model$markers[[node]]$gene
    
    x2 = x2/matrixStats::rowMaxs(x2, na.rm = T)
    
    pred.prob = predict(model$models[[node]], x2)
    pred.class = ifelse(pred.prob > 0.5, max(new_ids), min(new_ids))
    
    pred_df = data.frame(barcode = colnames(temp_seurat),
                         pred = pred.class,
                         prob = pred.prob,
                         census_clusters = temp_seurat$census_clusters,
                         u1 = temp_seurat@reductions$umap@cell.embeddings[,1],
                         u2 = temp_seurat@reductions$umap@cell.embeddings[,2],
                         stringsAsFactors = F)
    
    # if(!is.null(allowed_nodes) & any(pred_df$pred %in% allowed_nodes == F)){
    #   pred_df$pred = new_ids[new_ids %in% allowed_nodes]
    # }
    
    pred_res = contour_adjust(pred_df)
    
    return(pred_res)
  }
}

assignInNamespace("predict_node", predict_node_M, ns = "Census")

# Seurat
sample_id = 'I007_072'
obj = readRDS(paste0("./integrative_analysis/processed_obj/",sample_id,"_seurat_obj.rds"))

DefaultAssay(obj) = 'RNA'
obj = obj %>% NormalizeData() %>% ScaleData() %>%
  FindVariableFeatures() %>% RunPCA() %>% RunUMAP(dims=1:50)

# Census
res = census_main(obj, organ = 'Muscle')




# Plot Census annotations
p_census = ggplot(res$pred, aes(umap1, umap2, color = celltype)) + 
  geom_point(size = 0.7, shape = 16) + 
  geom_label_repel(data = res$pred %>% 
                     group_by(celltype) %>% 
                     summarize(u1 = mean(umap1), u2 = mean(umap2)), 
                   force = 100, 
                   aes(u1,u2, label=celltype), 
                   size = 4, color = 'black', 
                   label.padding = 0.1, 
                   max.overlaps = 20) +
  scale_color_manual(values=colorRampPalette(brewer.pal(10, "Spectral"))(length(unique(res$pred$celltype)))) +
  theme_classic() + 
  labs(title = paste0(sample_id,'- Census prediction'))+
  theme(legend.position = 'none', 
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(hjust = 0.5,face = 'bold'),
        axis.title = element_text(color = 'black', size = 7))

ggsave(paste0("./integrative_analysis/figures/",sample_id,'census_cell_type.pdf'),p_census,
       width = 6,height = 4)
