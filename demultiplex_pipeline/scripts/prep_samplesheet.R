
.libPaths('/home/i439h/projects/heroes-aya-pools/temp_analysis/tools/R/4.0/')

setwd('/home/i439h/projects/heroes-aya-pools/AG_Thongjuea/Code/10x_multiomics/Sadeghi/demultiplex_pipeline/test/')

#load libraries
library(stringr)
library(readr)
library(tidyr)

## we are going to have samplesheet dataframe with  columns PoolNumber, Sample, BAMPath, VCFPath, BarcodePath
#samplesheet = data.frame(PoolNumber= NA, Sample=NA, BAMPath=NA, VCFPath=NA, BarcodePath=NA)

## find vcf files
dir = '/home/i439h/projects/heroes-aya/AG_Thongjuea/Sub_project1/Result/bulkWGS/BJ-WGS'
files <- list.files(path =dir,pattern = 'snpEff.ann.vcf',full.names = T,recursive = T)
files = files[grep('control',files,invert = T)]

## extract sample names
pattern <- "I[0-9]+_[0-9]+"
samples = regmatches(files, regexpr(pattern, files))

samplesheet = data.frame(VCFPath = files, Sample = samples)

## PoolNumber, BamPath and BarcodePath must be added manually
## include only pool3 and 4 for test
# > AS-1058560-LR-69403 = I010_013
# > AS-1058563-LR-69404 = I005_014
# > AS-1058556-LR-69402 = I023_037
# > AS-1058553-LR-69401 = I020_001
# >
#   > pool4
# > AS-1058561-LR-69403 = I010_022
# > AS-1058553-LR-69401 = I020_004
# > AS-1058564-LR-69404 = I005_016
# > AS-1058556-LR-69402 = I023_047

pool3 = c('I010_013','I005_014','I023_037','I020_001')
pool4 = c('I010_022','I020_004','I005_016','I023_047')

samplesheet = subset(samplesheet,Sample %in% c(pool3,pool4))
samplesheet$PoolNumber = ifelse(samplesheet$Sample %in% pool3,'pool3','pool4')


## add Bam files paths
dir = "/home/i439h/projects/heroes-aya-pools/AG_Thongjuea/Result/10x_multiomics/CellRanger_ARC"
bamfiles = list.files(path =dir,pattern = 'gex_possorted_bam.bam',full.names = T,recursive = T)
bamfiles = bamfiles[grep('bai',bamfiles,invert = T)]
pool3Bam = bamfiles[grep('pool3',bamfiles)]
pool4Bam = bamfiles[grep('pool4',bamfiles)]

samplesheet$BAMPath = ifelse(samplesheet$PoolNumber == 'pool3',pool3Bam,pool4Bam)

## add barcode files paths
barcode = list.files(path =dir,pattern = 'barcodes.tsv.gz',full.names = T,recursive = T)
barcode = barcode[grep('filtered_feature_bc_matrix',barcode)]
pool3barcode = barcode[grep('pool3',barcode)]
pool4barcode = barcode[grep('pool4',barcode)]

samplesheet$BarcodePath = ifelse(samplesheet$PoolNumber == 'pool3',pool3barcode, pool4barcode)
samplesheet = samplesheet %>% dplyr::select( PoolNumber, Sample, BAMPath, VCFPath, BarcodePath)

## save file
write_tsv(samplesheet, file = 'samplesheet.tsv',quote = 'none')
