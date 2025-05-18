setwd("/omics/odcf/analysis/OE0290_projects/heroes-aya_pools/AG_Thongjuea/Code/10x_multiomics/Sadeghi/cellRanger_arc/")

pools= paste0('pool',9:33)
fastq = data.frame(fastqs =paste0(c('/omics/odcf/analysis/OE0290_projects/heroes-aya_pools/AG_Thongjuea/Dataset/10x_multiomics/linked_fastqs/scrna/',
                             '/omics/odcf/analysis/OE0290_projects/heroes-aya_pools/AG_Thongjuea/Dataset/10x_multiomics/linked_fastqs/atac/'),
                             rep(pools,1,each=2)),
                   sample = rep(pools,1,each=2),
                   library_type = c("Gene Expression","Chromatin Accessibility"))

for(p in pools){
  dat = subset(fastq,sample==p)
  write.csv(dat,file = paste0(p,'.csv'),quote = F,row.names = F)
}
