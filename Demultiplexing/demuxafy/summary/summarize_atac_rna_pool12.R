
cols = c("barcode",'assignment')
v_rna = read.delim(file ='./Demultiplexing_results/demuxafy/RNA/genotype/vireo/pool12/donor_ids.tsv',
                 header = T)
colnames(v_rna)[1:2] = cols
v_rna = v_rna[,cols]
colnames(v_rna) = paste0('vireo_rna_',colnames(v_rna))

v_atac = read.delim(file ='./Demultiplexing_results/demuxafy/ATAC/no_genotype/vireo/pool12/donor_ids.tsv',
                   header = T)
colnames(v_atac)[1:2] = cols
v_atac = v_atac[,cols]
colnames(v_atac) = paste0('vireo_atac_',colnames(v_atac))

soup_rna = read.delim(file ='./Demultiplexing_results/demuxafy/RNA/no_genotype/souporcell/pool12/clusters.tsv',
                      header = T)
soup_rna$assignment <- case_when(
  soup_rna$status != 'singlet' ~ soup_rna$status,
  soup_rna$status == 'singlet' ~ paste0('donor', soup_rna$assignment)
)

soup_rna = soup_rna[,cols]
colnames(soup_rna) = paste0('soup_rna_',colnames(soup_rna))

soup_atac = read.delim(file ='./Demultiplexing_results/demuxafy/ATAC/no_genotype/souporcell/pool12/clusters.tsv',
                      header = T)
soup_atac$assignment <- case_when(
  soup_atac$status != 'singlet' ~ soup_atac$status,
  soup_atac$status == 'singlet' ~ paste0('donor', soup_atac$assignment)
)

soup_atac = soup_atac[,cols]
colnames(soup_atac) = paste0('soup_atac_',colnames(soup_atac))

merged = cbind(v_rna,v_atac,soup_atac,soup_rna)

merged$vireo_atac_assignment = gsub('donor2',"I070_032", merged$vireo_atac_assignment)
merged$soup_atac_assignment = gsub('donor3',"I070_032", merged$soup_atac_assignment)
merged$soup_rna_assignment = gsub('donor0',"I070_032", merged$soup_rna_assignment)

merged_pass = merged %>% dplyr::select(vireo_rna_barcode,contains('assignment'))

# Add the majority vote column
merged_pass$majority_vote <- apply(merged_pass[, 2:5], 1, function(x) {
  # Get the most frequent value (excluding 'unassigned' or ties)
  votes <- table(x)
  votes <- votes[!names(votes) %in% "unassigned"]  # Exclude 'unassigned' from the vote count
  if(length(votes) == 0) return("unassigned")  # If all values are 'unassigned'
  majority <- names(votes)[votes == max(votes)]
  if(length(majority) > 1) return("tie") else return(majority) # Handle ties
})

# Display the updated data frame
print(merged_pass)
colnames(merged_pass)[1] = 'barcode'

write.table(merged_pass,'./Demultiplexing_results/demuxafy/summary/RNA_ATAC_pool12_majority.tsv',sep = '\t')
