
# vireoSummary.R
setwd("/omics/odcf/analysis/OE0536_projects/heroes-aya/scRNA-Seq/Demultiplexing_results/")
## load packages
pacman::p_load(ggplot2,Seurat, patchwork,ggthemes,GGally,ggpubr, Hmisc)

pools = paste0('pool',1:4)

## read vireo donor_id results
files = list.files(path = './vireo/results',pattern = "donor_id",full.names = T,recursive = T)

plots <- list()
for (i in 1:4) {
  dat = read.delim(file = files[i])
  
  p = ggplot(dat, aes(x = prob_max, y = prob_doublet, fill = donor_id))+
    geom_point(size=1, shape= 21,color = 'grey')+
    scale_fill_brewer(palette = 'Set1')+
    labs(x = "singlet prob",
         y = "doublet prob", 
         title = pools[i])+
    ggthemes::theme_clean()
  plots[[i]] = p
}

# merge plots
pp = plots[[1]]+plots[[2]]+plots[[3]]+plots[[4]]

# save plots
ggsave('./downstream_analysis/results/figures/vireo_pools_summary.pdf',plot = pp, width = 14,height = 8)

## summary stats
sums = list.files(path = './vireo/results',pattern = "summary",full.names = T,recursive = T)
dat = data.frame()
plt = list()

for(i in 1:4){
  input= read.delim(file = sums[i])
  input$pool = pools[i]
  dat = rbind(dat,input)
  
}
colnames(dat) = c("donor_id","cells",'pool')
## plot histogram 

p = ggplot(dat, aes(y = donor_id, x = cells, fill = donor_id))+
  geom_bar(stat = "identity",position = 'dodge',color = 'black',width = 0.7)+
  facet_wrap(~pool,scales = 'free')+
  geom_text(aes(label = cells),nudge_x =-1)+
  scale_x_log10()+
 #scale_fill_brewer(palette = 'Set1')+
  ggthemes::theme_clean(base_size = 12)

ggsave('./downstream_analysis/results/figures/vireo_stats.pdf',plot = p,width = 18,height = 7)
