
# souporcellSummary_atac.R
setwd("/home/i439h/projects/heroes-aya-pools/temp_analysis/scRNA-Seq/Demultiplexing_results/downstream_analysis/multi_omix/")
.libPaths("/home/i439h/projects/heroes-aya-pools/temp_analysis/tools/R/4.0/")

## load packages
pacman::p_load(ggplot2, patchwork,ggthemes,GGally, readr)

pools = paste0('pool',3:4)

## read souporcell cluster results
dir = "/home/i439h/projects/heroes-aya-pools/temp_analysis/scRNA-Seq/Demultiplexing_results/souporcell/results/multiomix/RNA/with_bulkWGS_BJ/"
files = list.files(path = dir, pattern = "clusters.tsv",full.names = T,recursive = T)
#files = files[grep("pool[0-9]_result",files)]

plots <- list()
allstat = onlysing = data.frame()
for (i in seq_along(pools)) {
  dat = read.delim(file = files[i])
  
  p = ggplot(dat, aes(x = log_prob_singleton, y = log_prob_doublet, fill = status))+
    geom_point(size=1, shape= 21,color = 'grey')+
    scale_fill_brewer(palette = 'Set1')+
    labs(x = " log singlet prob",
         y = "log doublet prob", 
         title = pools[i])+
    theme_clean()
  plots[[i]] = p
  
  ## tables
  stat= as.data.frame(table(dat$status))
  colnames(stat) = c("status",'cells')
  stat$pool = pools[i]
  allstat = rbind(allstat,stat)
  
  ## clusters freq in singlet
  singlet = subset(dat,status=='singlet')
  singstat = as.data.frame(table(singlet$assignment))
  colnames(singstat) = c("donor",'cells')
  singstat$pool = pools[i]
  onlysing = rbind(onlysing,singstat)
  
}

onlysing$donor = paste0('donor',onlysing$donor)
# merge plots
gg <- plots[[1]]+plots[[2]]

# save plots
ggsave('./results/figures/souporcell_summary_multiomix_RNA_WGS_BJ.pdf',plot = gg,
       width = 10,height = 5)

## bar chart
p1 = ggplot(allstat, aes(x = cells, y= status, fill = status))+
  geom_bar(stat = "identity",position = 'dodge',color = 'black',width = 0.7)+
  facet_wrap(~pool,scales = 'free')+
  geom_text(aes(label = cells),nudge_x =-1)+
  scale_x_log10()+
  scale_fill_brewer(palette = 'Paired')+
  ggthemes::theme_clean(base_size = 14)

p2 = ggplot(onlysing, aes(x = cells, y= donor, fill = donor))+
  geom_bar(stat = "identity",position = 'dodge',color = 'black',width = 0.7)+
  facet_wrap(~pool,scales = 'free')+
  geom_text(aes(label = cells),nudge_x =-1)+
  scale_x_log10()+
  scale_fill_brewer(palette = 'Set1')+
  ggthemes::theme_clean(base_size = 14)

gg= p1 / p2

ggsave('./results/figures/souporcell_stats_multiomix_RNA_WGS_BJ.pdf',plot = gg,
       width = 10,height = 5)
