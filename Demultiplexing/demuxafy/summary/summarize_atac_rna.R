# Set working directory and library paths
setwd("/home/i439h/projects/pools/AG_Thongjuea/Result/10x_multiomics/Demultiplexing_results/demuxafy/")
.libPaths("/omics/odcf/analysis/OE0290_projects/heroes-aya_pools/AG_Thongjuea/Software/10x_multiomics/R")

# Load packages
library(Seurat)
library(dplyr)

## read demultiplexing summaries

# Set the working directory
path <- getwd()

# Define vectors for data types, SNP, and tools
data_type <- c("RNA", "ATAC")
snp <- c("genotype", "no_genotype")

completeDat = data.frame()
# Loop through each data type
for(dat in data_type) {
  dir <- paste0(path, '/', dat)
  
  # Define the patterns to search for
  patterns <- c('summary.tsv', 'soupocell_summary.tsv')
  
  # List files matching any of the patterns
  names <- list.files(path = dir, pattern = "summary.tsv",recursive = TRUE)
  names = names[grep('demuxalot',names,invert = T)]
  files <- list.files(path = dir, pattern ='summary.tsv',full.names = T,recursive = TRUE)
  files = files[grep('demuxalot',files,invert = T)]
  
  # Print or process the files
  print(files)
  
  df <- data.frame(files = files, names = names, stringsAsFactors = FALSE)
  
  # Use the separate function to split the string
  df_separated <- df %>%
    separate(names, into = c("genotype", "tool", "pool", "file"), sep = "/", remove = TRUE)
  
  # Drop the file column as it's not needed
  df_separated <- df_separated %>%
    select(-file)
  df_separated$data_type = dat
  
  completeDat=rbind(completeDat,df_separated)
  
  # Prbind()# Print the result
  print(df_separated)
}

##genotype
dat_gene = subset(completeDat,genotype == 'genotype')
dat_gene = dat_gene[grep('combined',dat_gene$files,invert = T),]

alldat = data.frame()

  for( i in 1:nrow(dat_gene)){
    dat = read.delim(dat_gene$files[i])
    colnames(dat)= c('donor','assignment')
    dat$pool = dat_gene$pool[i]
    dat$tool = dat_gene$tool[i]
    dat$data_type = dat_gene$data_type[i]
    alldat = rbind(alldat,dat)
}

alldat <- alldat %>%
  mutate(donor = case_when(
    pool == 'pool3' & donor == "0" ~ "I005_014",
    pool == 'pool3' & donor == "1" ~ "I010_013",
    pool == 'pool3' & donor == "2" ~ "I020_001",
    pool == 'pool3' & donor == "3" ~ "I023_037",
    pool == 'pool4' & donor == "0" ~ "I005_016",
    pool == 'pool4' & donor == "1" ~ "I010_022",
    pool == 'pool4' & donor == "2" ~ "I020_004",
    pool == 'pool4' & donor == "3" ~ "I023_047",
    TRUE ~ donor # keep original donor if no match
  ))

#plots
ggplot(alldat,aes(x = donor,y = assignment, fill = data_type))+
  geom_bar(stat = 'identity',color = 'black',position = 'dodge')+
  facet_wrap(~pool,ncol = 1,scales = 'free')+
  scale_fill_brewer(palette = 'Accent')+
  ggthemes::theme_clean()


### no_genotype
completeDat = data.frame()


# Loop through each data type
for(dat in data_type) {
  dir <- paste0(path, '/', dat)
  
  # Define the patterns to search for
  patterns <- c('donor_ids.tsv', 'clusters.tsv')
  
  # List files matching any of the patterns
  names <- list.files(path = dir, pattern = paste0(patterns,collapse = '|'),recursive = TRUE)
  names = names[grep('demuxalot',names,invert = T)]
  files <- list.files(path = dir, pattern = paste0(patterns,collapse = '|'),full.names = T,recursive = TRUE)
  files = files[grep('demuxalot',files,invert = T)]
  
  # Print or process the files
  print(files)
  
  df <- data.frame(files = files, names = names, stringsAsFactors = FALSE)
  
  # Use the separate function to split the string
  df_separated <- df %>%
    separate(names, into = c("genotype", "tool", "pool", "file"), sep = "/", remove = TRUE)
  
  # Drop the file column as it's not needed
  df_separated <- df_separated %>%
    select(-file)
  df_separated$data_type = dat
  
  completeDat=rbind(completeDat,df_separated)
  completeDat = completeDat[grep('no_genotype',completeDat$files),]
  
  # Prbind()# Print the result
  print(df_separated)
}


## keep only pools with results
pools = paste0('pool',5:8)
completeDat = completeDat[completeDat$pool %in% pools, ]

dat_nogene = subset(completeDat,genotype == 'no_genotype')
dat_nogene = dat_nogene[grep('combined',dat_nogene$files,invert = T),]

alldat = data.frame()

for( i in 1:nrow(dat_nogene)){
  dat = read.delim(dat_nogene$files[i])
  if(dat_nogene$tool[i] == 'souporcell'){
  
  dat = dat %>% mutate(assignment = case_when(status == 'unassigned' ~ 'unassigned',
                               status == 'doublet'~ 'doublet',
                               assignment =='0'~'donor0',
                               assignment=='1'~ 'donor1',
                               assignment=='2'~ 'donor2',
                               assignment=='3'~ 'donor3',
                               TRUE ~ assignment))
  dat = dat[,c(1,3)]
  }else{
   dat = dat[,c(1,2)]  
  }
  colnames(dat)= c('donor','assignment')
  dat$pool = dat_nogene$pool[i]
  dat$tool = dat_nogene$tool[i]
  dat$data_type = dat_nogene$data_type[i]
  alldat = rbind(alldat,dat)
}

summary_data <- alldat %>%
  group_by(pool, tool, data_type, assignment) %>%
  summarise(count = n(), .groups = 'drop')

# Bar plot of assignments per pool, tool, and data type
p = ggplot(summary_data, aes(x = assignment, y = count, fill = data_type)) +
  geom_bar(stat = "identity",color ='black', position = position_dodge()) +
  facet_wrap(pool ~ tool, nrow = 2) +
  scale_fill_brewer(palette = 'Set1')+
  theme_minimal() +
  labs(title = "Number of Assignments per Pool, Tool, and Data Type",
       x = "Assignment",
       y = "Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save the plot
ggsave(filename = "./summary/assignments_per_pool_tool_data_type.pdf",width = 8,height = 5.5)

## stacked bar by tool
p2 = ggplot(summary_data, aes(x = data_type, y = count, fill = assignment)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap( pool ~ tool , scales = "free") +
  scale_fill_brewer(palette = 'Set1')+
  theme_minimal() +
  labs(title = "Stacked Number of Assignments by Tool",
       x = "Tool",
       y = "Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
# Save the plot
ggsave(filename = "./summary/stacked_assignments_per_pool_tool_data_type.pdf",width = 8,height = 5.5)

#pie chart

# Function to plot an enhanced pie chart for a specific pool
plot_pie_chart <- function(pool_name) {
  pool_data <- filter(summary_data, pool == pool_name)
  pool_data <- pool_data %>%
    group_by(tool, data_type, assignment) %>%
    summarise(count = sum(count), .groups = 'drop') %>%
    mutate(percentage = count / sum(count) * 100,
           label = paste0(assignment, "\n", round(percentage, 1), "%"))
  
  ggplot(pool_data, aes(x = "", y = count, fill = assignment)) +
    geom_bar(stat = "identity", width = 1, color = "white") +
    coord_polar("y") +
    facet_wrap(~ tool + data_type, scales = "free") +
    theme_minimal() +
    labs(title = paste("Distribution of Assignments in", pool_name),
         x = NULL, y = NULL) +
    theme(axis.text.x = element_blank(), axis.ticks = element_blank(),
          panel.grid = element_blank(), 
          plot.title = element_text(hjust = 0.5),
          strip.text = element_text(size = 8)) +
    geom_text_repel(aes(label = label), position = position_stack(vjust = 0.5), size = 3, show.legend = FALSE) +
    scale_fill_brewer(palette = 'Set1')+
    guides(fill = guide_legend(title = "Assignment Type"))
}


# Generate pie charts for each pool
pool_names <- unique(summary_data$pool)
for (pool in pool_names) {
  print(plot_pie_chart(pool))
  ggsave(paste0("./summary/pie_chart_", pool, ".pdf"),width = 10,height = 8)
}

