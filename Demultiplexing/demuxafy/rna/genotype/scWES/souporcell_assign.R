setwd("/home/i439h/projects/heroes-aya-pools/AG_Thongjuea/Result/10x_multiomics/Demultiplexing_results/demuxafy/RNA/genotype/scWES")
.libPaths(c('/home/i439h/projects/heroes-aya-pools/AG_Thongjuea/Software/10x_multiomics/R/'))
library(stringr)
library(readr)

# Diagnostic print to check received arguments
args <- commandArgs(trailingOnly = TRUE)
print(paste("Received args:", paste(args, collapse=", ")))

pool <- args[1]

print(paste("Pool:", pool))

# Iterate over specified pool only
clust_file <- paste0("./souporcell/", pool, "/clusters.tsv")
ids_file <- paste0("./souporcell/", pool, "/Genotype_ID_key.txt")

# Read clusters file
clust <- read_delim(clust_file)

# Read matched IDs
ids <- read_delim(ids_file)
ids$Cluster_ID <- 0:3
ids$Cluster_ID <- as.character(ids$Cluster_ID)
head(ids)

# Replace clusters with IDs
clust$assignment <- ids$Genotype_ID[match(clust$assignment, ids$Cluster_ID)]
clust$assignment[clust$status == 'doublet'] <- 'doublet'
clust$assignment[clust$status == 'unassigned'] <- 'unassigned'

# Write modified clusters file
output_file <- paste0("./souporcell/", pool, "/clusters_final.tsv")
write_tsv(clust, file = output_file)
