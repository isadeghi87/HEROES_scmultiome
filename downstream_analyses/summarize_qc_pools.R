## here we summarize the qc metrics for pooled sc-multiomics data 

setwd(dir = '/home/i439h/projects/heroes-aya-pools/AG_Thongjuea/Result/10x_multiomics/CellRanger_ARC/')

#libraries
library(stringr)

sumFiles = list.files(path = getwd(),pattern = "summary.csv",full.names = T,recursive = T)
sumFiles = sumFiles[grep('outs',sumFiles)]

pools = paste0('pool',1:8)

## iterate over pools

alldata = data.frame()
for (p in pools){
  csv = sumFiles[grep(p,sumFiles)]
  dat = read.csv(csv)
  dat$pool = p
  alldata = rbind(alldata,dat)
}


alldata$average.cells.per.sample = alldata$Estimated.number.of.cells/4
alldata = alldata[,-c(1:3)]

write.csv(alldata,"/home/i439h/projects/heroes-aya-pools/AG_Thongjuea/Result/10x_multiomics/summary/pools_qc_summary.csv")
