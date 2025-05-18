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
  pattern='combined_results_w_combined_assignments.tsv'
  
  # List files matching any of the patterns
  names <- list.files(path = path, pattern = pattern,recursive = TRUE)
  files <- list.files(path = path, pattern = pattern,full.names = T,recursive = TRUE)
  
  # Print or process the files
  df <- data.frame(files = files, names = names, stringsAsFactors = FALSE)
  
  # Use the separate function to split the string
  df_separated <- df %>%
    separate(names, into = c('data_type',"genotype", "combined", "pool", "file"), sep = "/", remove = TRUE)

pools = c('pool3','pool4')
gen = "genotype"

subdat = subset(df_separated,genotype ==gen)
dat.list = list()

for( p in pools){
pDat = subset(subdat,pool ==p)  
for( i in 1:nrow(pDat)){
  dat = read.delim(pDat$files[i])
  dat = dat[,c('Barcode','MajoritySinglet_Individual_Assignment')]
  data_type = pDat$data_type[i]
  colnames(dat) = paste0(data_type,'_',colnames(dat))
  dat$pool = p
 dat.list[[i]] <- dat
}

merged = do.call('cbind',dat.list)
merged = merged[,c(1,2,5,6)]
colnames(merged)[1] = 'barcode'

## final vote based on at least one assignment
merged$Donor_Assignment <- ifelse(
  merged$ATAC_MajoritySinglet_Individual_Assignment != "doublet",
  merged$ATAC_MajoritySinglet_Individual_Assignment,
  ifelse(
    merged$RNA_MajoritySinglet_Individual_Assignment != "doublet",
    merged$RNA_MajoritySinglet_Individual_Assignment,
    "doublet"
  )
)

write.table(merged,paste0('./summary/RNA_ATAC_',p,'_genotype_majority.tsv'),quote = F,row.names = F)
}

