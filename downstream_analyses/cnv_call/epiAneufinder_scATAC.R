setwd("/home/i439h/projects/")
.libPaths("/omics/odcf/analysis/OE0290_projects/heroes-aya_pools/AG_Thongjuea/Software/10x_multiomics/R/")

library(Rsamtools)
library(BSgenome)
library(epiAneufinder)

epiAneufinder(input="pools/AG_Thongjuea/Result/10x_multiomics/CellRanger_ARC/pool12/outs/individual_bams/", #Enter path to your fragments.tsv file or the folder containing bam files
              outdir="./pools/AG_Thongjuea/Result/10x_multiomics/cnv/results/integration/cnv/epiAneufinder_results", #Path to the directory where results should be written 
              blacklist="pools/AG_Thongjuea/Dataset/10x_multiomics/1k_genome/hg38-blacklist.v2.bed", #Path to bed file that contains the blacklisted regions of your genome
              windowSize=1e5, 
              genome="BSgenome.Hsapiens.UCSC.hg38", #Substitute with relevant BSgenome
              exclude=c('chrX','chrY','chrM'), 
              reuse.existing=TRUE,
              title_karyo="Karyogram of sample I070_032", 
              ncores=4,
              minFrags=20000,
              minsizeCNV=0,
              k=4,
              plotKaryo=TRUE)


#Load result table
res_table<-read.table("./pools/AG_Thongjuea/Result/10x_multiomics/cnv/results/integration/cnv/epiAneufinder_results/epiAneufinder_results/results_table.tsv")

#Need to reformat column names as dash are replaced by points in R
colnames(res_table) <-gsub("\\.","-",colnames(res_table))

res_table = res_table[,c(1:4,6)]

subclones<-split_subclones(res_table,tree_depth=2,
                           plot_tree=TRUE,
                           plot_path="pools/AG_Thongjuea/Result/10x_multiomics/cnv/results/integration/cnv/epiAneufinder_results/epiAneufinder_results/subclones.pdf",
                           plot_width=4,
                           plot_height=3)

#Data frame with the subclone split
head(subclones)
