setwd("/home/i439h/projects/heroes-aya-pools/temp_analysis/scRNA-Seq/Demultiplexing_results/downstream_analysis/")
.libPaths("/home/i439h/projects/heroes-aya-pools/temp_analysis/tools/R/4.0/")

## load packages
library(ggplot2); library(patchwork);library(ggthemes);library(GGally);library(readr)
library(VennDiagram);library(stringr)

pools = paste0('pool',3:4)
lib= c("vireo",'Souporcell')


## read vireo donor_id results
dir = "/home/i439h/projects/heroes-aya-pools/temp_analysis/scRNA-Seq/Demultiplexing_results/vireo/results/"
files = list.files(path = dir,pattern = "donor_id",full.names = T,recursive = T)
files = files[grep('WGS_BJ',files)]
vireo = data.frame()
for(p in seq_along(pools)){ 
  dat = read.delim(file = files[p])
  dat$pool=pools[p]
  vireo=rbind(vireo,dat)  
}
colnames(vireo) = str_replace_all(colnames(vireo),c("cell"="barcode",'donor_id'='assignment'))
vireo = subset(vireo, !assignment %in% c("doublet","unassigned"))



## souporcell
## read souporcell cluster results
dir = "/home/i439h/projects/heroes-aya-pools/temp_analysis/scRNA-Seq/Demultiplexing_results/souporcell/results/multiomix/RNA/with_bulkWGS_BJ/"
files = list.files(path = dir, pattern = "clusters.tsv",full.names = T,recursive = T)

soup = data.frame()
for (i in seq_along(pools)) {
  dat = read.delim(file = files[i])
  dat$pool=pools[i]
  soup=rbind(soup,dat)  
}
soup = subset(soup,status=='singlet')  
soup = soup[,c('barcode','assignment','pool')]
soup$lib= 'souporcell'
vireo = vireo[,c('barcode','assignment','pool')]
vireo$lib= 'vireo'
alldat = rbind(soup,vireo)


alldat$assignment  = str_replace_all(alldat$assignment,
                                     c('AS-1058563-LR-69404'="0",
                                       'AS-1058564-LR-69404'="0",
                                       'AS-1058560-LR-69403'="1",
                                       'AS-1058561-LR-69403'="1",
                                       'AS-1058553-LR-69401'="2",
                                       'AS-1058553-LR-69401'="2",
                                       'AS-1058556-LR-69402'="3",
                                       'AS-1058556-LR-69402'="2"))




library(RColorBrewer)
donor=unique(alldat$assignment)

for(p in pools){
  dat = subset(alldat,pool==p)
  for(d in donor){
    # filter each donor
    datDonor= subset(dat,assignment==d)
    
    #scRNA
    soup=subset(datDonor,lib=='souporcell')$barcode
    #multiomix
    vireo=subset(datDonor,lib=='vireo')$barcode
    int = length(intersect(soup,vireo))
    v = VennDiagram::venn.diagram(x = list(soup,vireo),filename = paste0("./multi_omix/results/figures/comparison/multiomics_",p,"_donor",d,"_venn.png"),
                                  output=FALSE,category.names = c("souporcell","Vireo"),
                                  fill=c("red",'blue'),
                                  main = paste0(p,'_donor',d),
                                  # Output features
                                  imagetype="png")
  }
}
