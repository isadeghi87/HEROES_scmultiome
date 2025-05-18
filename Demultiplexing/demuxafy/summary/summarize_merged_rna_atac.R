# Set working directory and library paths
## here we compare the demultiplexing results from merged (RNA+ATAC) vc RNA and ATAC
setwd("/home/i439h/projects/pools/AG_Thongjuea/Result/10x_multiomics/Demultiplexing_results/demuxafy/")
.libPaths("/omics/odcf/analysis/OE0290_projects/heroes-aya_pools/AG_Thongjuea/Software/10x_multiomics/R")

# Load packages
library(Seurat)
library(dplyr)

## read demultiplexing summaries

# Set the working directory
path <- getwd()

gen = c('genotype')
data_type = c('ATAC','RNA','merged')
pools = c('pool3','pool4')
alldat = amb_dat = data.frame()

for(p in pools){
  for (g in gen){
    for(d in data_type){
      dat = read.delim(file = paste0(d,"/",g,"/souporcell/",p,"/souporcell_summary.tsv"))
      dat$genotype = g
      dat$data = d
      dat$pool = p
      amb = read.table(file = paste0(d,"/",g,"/souporcell/",p,"/ambient_rna.txt"))
      ambient_RNA =  as.numeric(gsub("%", "", amb$V5))
      amb_dat = rbind(amb_dat,cbind(g,d,p,ambient_RNA))
      alldat = rbind(alldat,dat)
    }  
  }
}

colnames(alldat)[1:2] = c('donor','assignments')
alldat <- alldat %>%
  mutate(donor = case_when(
    donor ==0 ~ "donor0"  ,
    donor == 1 ~ "donor1" ,
    donor ==2 ~"donor2" ,
    donor ==3 ~ "donor3" ,
    TRUE ~ donor # keep original donor if no match
  ))


singlets = alldat
donors =paste0('donor',0:3)
singlets$donor[singlets$donor %in% donors] ='singlet'

# Calculate the sums for each assignment in each data type
result <- aggregate(assignments ~ donor + data+genotype+pool, data = singlets, sum)
result$donor =factor(result$donor,levels=c('singlet','doublet','unassigned'))

(pp = ggplot(result, aes(x = donor, y = assignments, fill = data,label=assignments)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap(~pool, scales = "free") +
    scale_fill_brewer(palette = 'Set1')+
    scale_y_log10()+
    theme_bw() +
    labs(title = "merged (RNA+ATAC) vs separate RNA and ATAC",subtitle = 'using WGS data') +
    geom_text(aes(label = assignments),position = position_dodge(width = 1),vjust =-0.5)+
    theme(axis.text.x = element_text(angle = 45, hjust = 1)))

ggsave(filename = "./summary/merged_vs_RNA_ATAC_barplot.pdf",plot = pp,width = 8,height = 5)

## ambient RNA plot 
amb_dat$ambient_RNA = as.numeric(amb_dat$ambient_RNA)

(amb_p = ggplot(amb_dat,aes(x=d,y=ambient_RNA,fill=d))+
    geom_bar(stat = 'identity')+
    facet_wrap(~p)+
    scale_fill_brewer(palette = 'Set1')+
    labs(x = ''))

ggsave(filename = "./summary/ambient_RNA_barplot.pdf",plot = amb_p, width = 8,height = 5)
