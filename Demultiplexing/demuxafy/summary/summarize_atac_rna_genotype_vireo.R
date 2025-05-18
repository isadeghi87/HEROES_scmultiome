# Set working directory and library paths
setwd("/home/i439h/projects/pools/AG_Thongjuea/Result/10x_multiomics/Demultiplexing_results/demuxafy/")
.libPaths("/omics/odcf/analysis/OE0290_projects/heroes-aya_pools/AG_Thongjuea/Software/10x_multiomics/R")

# Load packages
library(Seurat)
library(dplyr)
library(tidyr)

## read demultiplexing summaries

# Set the working directory
path <- getwd()

# Define vectors for data types, SNP, and tools
data_type <- c("RNA", "ATAC")

completeDat = data.frame()
  
  # Define the patterns to search for
  pattern='donor_ids.tsv'
  
  # List files matching any of the patterns
  names <- list.files(path = path, pattern = pattern,recursive = TRUE)
  files <- list.files(path = path, pattern = pattern,full.names = T,recursive = TRUE)
  
  # Print or process the files
  df <- data.frame(files = files, names = names, stringsAsFactors = FALSE)
  
  # Use the separate function to split the string
  df_separated <- df %>%
    separate(names, into = c('data_type',"genotype", "combined", "pool", "file"), sep = "/", remove = TRUE)

pools = c('pool14')
gen = "genotype"

subdat = subset(df_separated,genotype ==gen & pool ==pools)
alldat =list()

for( i in 1:nrow(subdat)){
  dat = read.delim(subdat$files[i])
  dat = dat[,c(1:2)]
  data_type = subdat$data_type[i]
  colnames(dat) = paste0(data_type,'_',colnames(dat))
  dat$pool = subdat$pool[i]
 alldat[[i]] = dat
}

merged = do.call('cbind',alldat)
merged = merged[,c(1,2,3,5)]
colnames(merged)[1] = 'barcode'

## final vote based on at least one assignment
merged$Donor_Assignment <- ifelse(
  merged$ATAC_donor_id != "doublet",
  merged$ATAC_donor_id,
  ifelse(
    merged$RNA_donor_id != "doublet",
    merged$RNA_donor_id,
    "doublet"
  )
)

write.table(merged,paste0('./summary/RNA_ATAC_',p,'_genotype_vireo.tsv'),quote = F,row.names = F)

