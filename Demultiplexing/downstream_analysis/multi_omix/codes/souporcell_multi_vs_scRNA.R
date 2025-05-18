setwd("/home/i439h/projects/heroes-aya-pools/temp_analysis/scRNA-Seq/Demultiplexing_results/downstream_analysis/")
.libPaths("/home/i439h/projects/heroes-aya-pools/temp_analysis/tools/R/4.0/")

## load packages
library(ggplot2); library(patchwork);library(ggthemes);library(GGally);library(readr)
library("ggVennDiagram")

pools = paste0('pool',3:4)
lib= c("scRNAseq",'multiomix')

allstat = onlysing = data.frame()
for(l in lib){
## read souporcell cluster results
rna_dir = paste0("/home/i439h/projects/heroes-aya-pools/temp_analysis/scRNA-Seq/Demultiplexing_results/souporcell/results/",l)
files = list.files(path = rna_dir, pattern = "clusters.tsv",full.names = T,recursive = T)
files = files[grep("WGS_BJ",files)]
for (i in seq_along(pools)) {
  dat = read.delim(file = files[i])
  dat$pool = pools[i]
  dat$lib = l
  allstat = rbind(allstat,dat)
  
  # ## clusters freq in singlet
  # singlet = subset(dat,status=='singlet')
  # singstat = as.data.frame(table(singlet$assignment))
  # colnames(singstat) = c("donor",'cells')
  # singstat$pool = pools[i]
  # onlysing = rbind(onlysing,singstat)
  
}
}

allstat = subset(allstat,status=="singlet")
donor=unique(allstat$assignment)

library(RColorBrewer)

for(p in pools){
  dat = subset(allstat,pool==p)
for(d in donor){
  # filter each donor
  datDonor= subset(dat,assignment==d)
  
  #scRNA
  rna=subset(datDonor,lib=='scRNAseq')$barcode
  #multiomix
  multi=subset(datDonor,lib=='multiomix')$barcode
  int = length(intersect(rna,multi))
  v = VennDiagram::venn.diagram(x = list(rna,multi),filename = paste0("./multi_omix/results/figures/",p,"_donor",d,"_venn.png"),
                                output=FALSE,category.names = c("scRNA","Multiomix"),
                                fill=c("orange",'blue'),
                                main = paste0(p,'_donor',d),
                                # Output features
                                imagetype="png")
}
}
