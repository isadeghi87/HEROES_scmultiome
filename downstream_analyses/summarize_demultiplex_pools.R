
setwd("/omics/odcf/analysis/OE0290_projects/heroes-aya_pools/AG_Thongjuea/Result/10x_multiomics/Demultiplexing_results/demuxafy/RNA/")
library(stringr);library(ggplot2)

## list demultiplex files for all pools

links = list.files(path = getwd(),pattern = 'souporcell_summary.tsv',recursive = T,full.names = T)

pools = paste0('pool',1:8)

alldat = sumdat = data.frame()

for (p in pools){
  
  dem_file = links[grep(p,links)]
  dat = read.delim(dem_file)
  dat$pool = p
  if (nrow(dat)==6) {
    dat$Classification[1:4] = paste0('sample',1:4)
  } else {
    dat$Classification[1:3] = paste0('sample',1:3)
  }
  alldat = rbind(alldat,dat)
  # ## sum of singlets
  # singlet = sum(dat$Assignment.N[1:4])
  # singletDat = data.frame(Singlets = singlet,
  #                         doublet = dat$Assignment.N[dat$Classification=='doublet'],
  #                         unassigned = dat$Assignment.N[dat$Classification=='unassigned'],
  #                         pool=p)
  # sumdat = rbind(sumdat,singletDat)
}

colnames(alldat) = c('sample','cells','pool')

p1 = ggplot(alldat, aes(x = pool,y =cells ,fill = sample))+
  geom_bar(stat = 'identity',position = 'dodge',color ='black')+
  scale_fill_brewer(palette = 'Set3')+
  ggtitle('Demultiplexing results')

p2 = ggplot(alldat, aes(x = sample,y =cells ,fill = pool))+
  geom_bar(stat = 'identity',position = 'dodge',color ='black')+
  scale_fill_brewer(palette = 'Set3')+
  ggtitle('Demultiplexing results')

ggsave("../summary/RNA_pools_barplot.pdf",p1,width = 10,height = 6)
ggsave("../summary/RNA_pools_barplot2.pdf",p2,width = 10,height = 6)
write.csv(x = alldat,file = "../summary/RNA_pools_summary.csv")
