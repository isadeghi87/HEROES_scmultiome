
# vireoSummary.R
# souporcellSummary_atac.R
setwd("/home/i439h/projects/heroes-aya-pools/temp_analysis/scRNA-Seq/Demultiplexing_results/downstream_analysis/multi_omix/")
.libPaths("/home/i439h/projects/heroes-aya-pools/temp_analysis/tools/R/4.0/")
pacman::p_load(ggplot2, patchwork,ggthemes,GGally,ggpubr, Hmisc)

pools = paste0('pool',3:4)

## read vireo donor_id results
dir = "/home/i439h/projects/heroes-aya-pools/temp_analysis/scRNA-Seq/Demultiplexing_results/vireo/results/"
files = list.files(path = dir,pattern = "donor_id",full.names = T,recursive = T)
files = files[grep('multiomics',files)]

plots <- list()
for (i in seq_along(pools)) {
  dat = read.delim(file = files[i])
  dat$donor_id = str_replace_all(dat$donor_id,
                                 c(I005_014='AS-1058563-LR-69404',
                                   I005_016='AS-1058564-LR-69404',
                                   I010_013='AS-1058560-LR-69403',
                                   I010_022='AS-1058561-LR-69403',
                                   I020_001='AS-1058553-LR-69401',
                                   I020_004='AS-1058553-LR-69401',
                                   I023_037='AS-1058556-LR-69402',
                                   I023_047='AS-1058556-LR-69402'))
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
pp = plots[[1]]+plots[[2]]

# save plots
ggsave('./results/figures/vireo_summary_multiomix_RNA_WGS_BJ.pdf',plot = pp, width = 14,height = 8)

## summary stats
sums = list.files(path = dir, pattern = "summary.tsv",full.names = T,recursive = T)
sums = sums[grep('WGS_BJ',sums)]
dat = data.frame()
plt = list()

for(i in seq_along(pools)){
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

ggsave('./results/figures/vireo_stats_multiomix_RNA_WGS_BJ.pdf',plot = p,width = 18,height = 7)

res = 
