setwd("/home/i439h/projects/heroes-aya-pools/AG_Thongjuea/Result/10x_multiomics/Demultiplexing_results/")
.libPaths("/home/i439h/projects/heroes-aya-pools/temp_analysis/tools/R/4.0/")

## load packages
library(ggplot2); library(patchwork);library(ggthemes);library(GGally);library(readr)
library(VennDiagram);library(stringr)
library(dplyr);library(purrr)
library(tidyr)
pools = paste0('pool',3:4)
lib= c("vireo",'souporcell')


  dir = paste0("/home/i439h/projects/heroes-aya-pools/temp_analysis/scRNA-Seq/Demultiplexing_results/vireo/results/")
  files = list.files(path = dir,pattern = "donor_id",full.names = T,recursive = T)
  files = files[grep('multiomix',files)]
  files = files[grep('RNA',files)]
  
  
  # Read each file and keep only first 2 columns
  dfs <- lapply(files, function(file) {
    df <- read_tsv(file)[,1:2]
    colnames(df) = c('cell','donor_id')
    df$source <- basename(dirname(file))  # Adding a column to keep track of the source (like unique_snps, q2000_dp100, etc.)
    method = stringr::str_split(file,pattern = "/",simplify = T)
    method = method[ncol(method)-2]
    df$method <- method
    df$tool = 'vireo'
    return(df)
  })
  
  # Get counts for each donor from each file
  counts_list <- lapply(dfs, function(df) {
    table(df$donor_id)
  })
  
  # View counts
  counts_list
  
  # Venn diagram for cell overlap (assuming only considering 3 files for simplicity; you can adjust)
  # venn.diagram(
  #   x = list(
  #     dfs[[1]]$cell,
  #     dfs[[2]]$cell,
  #     dfs[[3]]$cell
  #   ),
  #   category.names = c("pool1", "pool2", "pool3"),
  #   output = TRUE
  # )
  
  # Matrix for cell assignments
  vireo_df <- bind_rows(dfs)
  
  #vireo <- spread(vireo_df, key = method, value = donor_id)
  
  library(dplyr)
  
  # Define a data frame for replacements
  replacements <- data.frame(
    source = rep(c('pool4_1016717', 'pool3_1016715'), each = 4),
    old_id = c('AS-1058561-LR-69403', 'AS-1058553-LR-69401', 'AS-1058564-LR-69404', 'AS-1058556-LR-69402',
               'AS-1058560-LR-69403', 'AS-1058563-LR-69404', 'AS-1058556-LR-69402', 'AS-1058553-LR-69401'),
    new_id = rep(c('donor0', 'donor1', 'donor2', 'donor3'), 2)
  )
  
  # Perform the replacements
  vireo_df <- vireo_df %>%
    left_join(replacements, by = c("source", "donor_id" = "old_id")) %>%
    mutate(donor_id = ifelse(!is.na(new_id), new_id, donor_id)) %>%
    dplyr::select(-new_id)
  
  # vireo_df$donor_id[vireo_df$source=='pool3_1016715'& vireo_df$donor_id=='AS-1058560-LR-69403']= 'donor0'
  # vireo_df$donor_id[vireo_df$source=='pool3_1016715'& vireo_df$donor_id=='AS-1058563-LR-69404']= 'donor1'
  # vireo_df$donor_id[vireo_df$source=='pool3_1016715'& vireo_df$donor_id=='AS-1058556-LR-69402']= 'donor2'
  # vireo_df$donor_id[vireo_df$source=='pool3_1016715'& vireo_df$donor_id=='AS-1058553-LR-69401']= 'donor3'
  
  # vireo_df$donor_id[vireo_df$source=='pool4_1016717'& vireo_df$donor_id=='AS-1058561-LR-69403']= 'donor0'
  # vireo_df$donor_id[vireo_df$source=='pool4_1016717'& vireo_df$donor_id=='AS-1058553-LR-69401']= 'donor1'
  # vireo_df$donor_id[vireo_df$source=='pool4_1016717'& vireo_df$donor_id=='AS-1058564-LR-69404']= 'donor2'
  # vireo_df$donor_id[vireo_df$source=='pool4_1016717'& vireo_df$donor_id=='AS-1058556-LR-69402']= 'donor3'
  # 
  
  vireo = spread(vireo_df,key = method,value = donor_id)
  
  
  #### souporcell ####
  dir = paste0("/home/i439h/projects/heroes-aya-pools/temp_analysis/scRNA-Seq/Demultiplexing_results/souporcell/results/")
  files = list.files(path = dir,pattern = 'clusters.tsv',full.names = T,recursive = T)
  files = files[grep('multiomix',files)]
  files = files[grep('RNA',files)]
  
  
  # Read each file and keep only first 2 columns
  dfs_soup <- lapply(files, function(file) {
    df <- read_tsv(file)[,c(1:3)]
    colnames(df) = c('cell','status','donor_id')
    df$donor_id[df$status=='unassigned']='unassigned'
    df$donor_id[df$status=='doublet']='doublet'
    df = df[,c(1,3)]
    df$source <- basename(dirname(file))  # Adding a column to keep track of the source (like unique_snps, q2000_dp100, etc.)
    method = stringr::str_split(file,pattern = "/",simplify = T)
    method = method[ncol(method)-2]
    df$method <- method
    df$tool = tool
    return(df)
  })
  
  # Get counts for each donor from each file
  counts_list <- lapply(dfs_soup, function(df) {
    table(df$donor_id)
  })
  
  
  # Matrix for cell assignments
  souporcell_df <- bind_rows(dfs_soup)[-1,]
  souporcell_df <- souporcell_df %>%
    mutate(donor_id = ifelse(!donor_id %in% c('doublet', 'unassigned'), 
                             paste0('donor', donor_id), 
                             donor_id))

  
  # Create the mapping table
  mapping <- data.frame(
    source = c(rep("pool4_1016717", 4), rep("pool3_1016715", 4)),
    method = rep("no_genotype", 8),
    original_id = c('donor0', 'donor1', 'donor2', 'donor3', 'donor0', 'donor1', 'donor2', 'donor3'),
    new_id = c('donor3', 'donor2', 'donor0', 'donor1', 'donor2', 'donor1', 'donor3', 'donor0')
  )
  
  # Apply the mapping to your dataframe
  souporcell_df <- souporcell_df %>%
    left_join(mapping, by = c("source", "method", "donor_id" = "original_id")) %>%
    mutate(
      donor_id = if_else(!is.na(new_id), new_id, donor_id)
    ) %>%
    dplyr::select(-new_id)  # Remove the new_id column
  
  # souporcell_df <- spread(souporcell_df, key = method, value = donor_id)
  allres = rbind(vireo_df,souporcell_df)
  donor = unique(souporcell_df$donor_id)
  pools = unique(souporcell_df$source)
  methods = unique(souporcell_df$method)
  
      for(m in methods){
        dat = subset(souporcell_df,method==m)
        tab = table(dat$donor_id,dat$source)
      }

  #### venn diagram
  binary_df <- souporcell_df %>%
    mutate(
      merged = ifelse(!is.na(merged) & merged != "doublet", 1, 0),
      no_genotype = ifelse(!is.na(no_genotype) & no_genotype != "doublet", 1, 0),
      q2000_dp100 = ifelse(!is.na(q2000_dp100) & q2000_dp100 != "doublet", 1, 0),
      unique_snps = ifelse(!is.na(unique_snps) & unique_snps != "doublet", 1, 0)
    )

  
  library(VennDiagram)
  
  venn_data <- list(
    merged = binary_df$cell[binary_df$merged == 1],
    no_genotype = binary_df$cell[binary_df$no_genotype == 1],
    q2000_dp100 = binary_df$cell[binary_df$q2000_dp100 == 1],
    unique_snps = binary_df$cell[binary_df$unique_snps == 1]
  )
  
  venn.plot <- venn.diagram(
    x = venn_data,
    category.names = c("merged", "no_genotype", "q2000_dp100", "unique_snps"),
    output = TRUE, fill=c("red",'blue','green','orange'),
      filename = './downstream_analysis/multi_omix/results/figures/souporcell_all.png'
  )
  
  souporcell =  spread(souporcell_df, key = method, value = donor_id)
  counts = souporcell_df %>% group_by(source,donor_id,method) %>% summarise(n=n())
  alldat = rbind(vireo_df,souporcell_df)
  
  ## match donors between souporcell and vireo
  library(dplyr)
  
  # Define a data frame for replacements
  replacements <- data.frame(
    tool = rep('souporcell', 8),
    source = rep(c('pool4_1016717', 'pool3_1016715'), each = 4),
    old_id = rep(c('donor0', 'donor1', 'donor2', 'donor3'), 2),
    new_id = c('donor2', 'donor0', 'donor1', 'donor3', 'donor1', 'donor0', 'donor3', 'donor2')
  )
  
  # Perform the replacements
  alldat <- alldat %>%
    left_join(replacements, by = c("tool", "source", "donor_id" = "old_id")) %>%
    mutate(donor_id = ifelse(!is.na(new_id), new_id, donor_id)) %>%
    dplyr::select(-new_id)
  alldat$method = str_replace_all(alldat$method,
                  c('merged'='WGS_1M_snps',
                    'q2000_dp100'= 'WGS_500k_snps',
                    'unique_snps'='WGS_100k_unique_snps'))
  # transform the data to make a matrix
  alldat_tf =  spread(alldat, key = method, value = donor_id)
  alldat_pool3_4 = subset(alldat,source %in% c('pool3_1016715','pool4_1016717'))
  alldat_pool3_4$method = paste0(alldat_pool3_4$tool,'_',alldat_pool3_4$method)
  alldat_pool3_4 %<>% dplyr::select(cell,source,donor_id,method) 
  alldat_tf_pool3_4 =  spread(alldat_pool3_4, key = method, value = donor_id)
  write_tsv(alldat_tf_pool3_4,file = './downstream_analysis/multi_omix/results/demultiplex_summary_matrix_pool3_4.tsv',col_names = T)
  
  ### mapping back to sample IDs
  # pool3
  # > AS-1058560-LR-69403 = I010_013 = donor0
  # > AS-1058563-LR-69404 = I005_014= donor1
  # > AS-1058556-LR-69402 = I023_037= donor2
  # > AS-1058553-LR-69401 = I020_001= donor3
  # >
  #   > pool4
  # > AS-1058561-LR-69403 = I010_022 = donor0
  # > AS-1058553-LR-69401 = I020_004 = donor1
  # > AS-1058564-LR-69404 = I005_016 = donor2
  # > AS-1058556-LR-69402 = I023_047= donor3
  
  replacements <- data.frame(
    source = rep(c('pool3_1016715','pool4_1016717'), each = 4),
    old_id = rep(c('donor0', 'donor1', 'donor2', 'donor3'), 2),
    new_id = c('I010_013', 'I005_014', 'I023_037', 'I020_001',
               'I010_022', 'I020_004', 'I005_016', 'I023_047'))
  
  # Perform the replacements
  alldat_ID <- alldat_pool3_4 %>%
    left_join(replacements, by = c("source", "donor_id" = "old_id")) %>%
    mutate(donor_id = ifelse(!is.na(new_id), new_id, donor_id)) %>%
    dplyr::select(-new_id)
  
  alldat_ID_tf = spread(alldat_ID,key = method, value = donor_id)
  write_tsv(alldat_ID_tf,file = './downstream_analysis/multi_omix/results/demultiplex_summary_matrix_pool3_4_ID.tsv',col_names = T)
  
  ## save one with 
  wgdata = alldat_ID_tf %>% dplyr::select(cell,source,vireo_WGS_1M_snps)
  write_tsv(wgdata,file = './downstream_analysis/multi_omix/results/demultiplex_vireo_wgs.tsv',col_names = T)
  
  
  
  ## compute the majority vote
  find_most_frequent <- function(df) {
    df %>%
      rowwise() %>%
      mutate(majority = {
        words <- c_across(-cell)  # Adjust column selection as needed, excluding 'cell'
        freqs <- table(words)
        most_freq <- names(freqs)[which.max(freqs)]
        if (length(most_freq) > 1) most_freq <- most_freq[1]  # Handle ties
        most_freq
      }) %>%
      ungroup()
  }
  
  # Example usage
   new_df <- find_most_frequent(alldat_ID_tf)
   write_tsv(new_df,'./downstream_analysis/multi_omix/results/demultiplex_majority_vote.tsv',col_names = T)
   
   pools = c('pool3_1016715','pool4_1016717')
  for (p in pools){
    dat = subset(new_df, source ==p)
    dat = dat %>% dplyr::select(cell,majority)
    write_tsv(dat,file = paste0('./downstream_analysis/multi_omix/results/',p,'_majority.tsv'),col_names = T)
  }
   
   
  write_tsv(allcounts,file = './downstream_analysis/multi_omix/results/demultiplex_allcounts.tsv',col_names = T)
  
  ## filter pool3&4
  subcount = subset(result,source %in% c('pool3_1016715','pool4_1016717'))
  
  # plot
  dodge_width <- 0.9
  pl = ggplot(subcount, aes(x = donor_id, y = n)) +
    geom_bar(aes(fill = method),width = 0.8, color = 'black', stat = 'identity', position = position_dodge(width = dodge_width)) +
    facet_grid(tool~source, scales = 'free') +
    scale_fill_brewer(palette = 'Set1') +
    geom_text(aes(label = n, group = method), size = 2, vjust = -0.25, position = position_dodge(width = dodge_width)) +
    theme_bw()
pl
  ggsave(plot = pl,filename = './downstream_analysis/multi_omix/results/figures/overall_summary_barplot.pdf',
         width = 14,height = 8)  
  
pool4 = subset(alldat_ID_tf,source=='pool4_1016717')  
library(knitr)
library(kableExtra)
library(webshot)

colorful_table <- kable(head(pool4), format = "html") %>%
  kable_styling(bootstrap_options = "striped", full_width = FALSE) %>%
  row_spec(0, background = "lightblue") %>%
  column_spec(1, background = "pink")



