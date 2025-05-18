library(stringr)

# Define the directory containing the original FASTQ files
original_dir <- "/omics/odcf/project/OE0290/heroes-aya_pools/sequencing/"

# List all FASTQ files recursively
all.f <- list.files(original_dir, 'fastq.gz', full.names = TRUE, all.files = TRUE, recursive = TRUE)

# Filter out files that do not contain 'core'
all.f <- all.f[grep('core', all.f, invert = TRUE)]

# Extract file paths containing 'sequence'
file_paths <- all.f[grep('sequence', all.f)]

# Define pools to include
pools <- paste0('pool', 11:33)

# Create a pattern to filter the pools
pattern <- paste0(pools, collapse = '|')

# Filter file paths based on the pool pattern
file_paths <- file_paths[grep(pattern, file_paths)]

# Function to create symbolic links
create_symlink_with_rename <- function(original_file_paths, link_base_dir) {
  # Initialize a list to store unique pool and data type combinations
  pool_data_type_map <- list()
  
  # Group files by pool and data type
  for (original_path in original_file_paths) {
    pool_number_match <- regmatches(original_path, regexpr("pool\\d+", original_path))
    pool_number <- if (length(pool_number_match) > 0) pool_number_match[1] else NA
    
    data_type_match <- regmatches(original_path, regexpr("10x_multiome_\\w+_sequencing", original_path))
    data_type <- if (length(data_type_match) > 0) sub("10x_multiome_", "", sub("_sequencing", "", data_type_match[1])) else NA
    
    sample_number_match <- regmatches(original_path, regexpr("AS-\\d+", original_path))
    sample_number <- if (length(sample_number_match) > 0) as.integer(sub("AS-", "", sample_number_match[1])) else NA
    
    lane_number_match <- regmatches(original_path, regexpr("LR-\\d+", original_path))
    lane_number <- if (length(lane_number_match) > 0) as.integer(sub("LR-", "", lane_number_match[1])) else NA
    
    if (!is.na(pool_number) && !is.na(data_type) && !is.na(sample_number) && !is.na(lane_number)) {
      pool_data_type <- paste(pool_number, data_type, sep = "_")
      if (is.null(pool_data_type_map[[pool_data_type]])) {
        pool_data_type_map[[pool_data_type]] <- list(files = c(), sample_numbers = c(), lane_numbers = c(), sample_map = list(), lane_map = list())
      }
      pool_data_type_map[[pool_data_type]]$files <- c(pool_data_type_map[[pool_data_type]]$files, original_path)
      pool_data_type_map[[pool_data_type]]$sample_numbers <- c(pool_data_type_map[[pool_data_type]]$sample_numbers, sample_number)
      pool_data_type_map[[pool_data_type]]$lane_numbers <- c(pool_data_type_map[[pool_data_type]]$lane_numbers, lane_number)
    }
  }
  
  # Assign sample and lane numbers for each pool and data type
  for (pool_data_type in names(pool_data_type_map)) {
    file_paths <- pool_data_type_map[[pool_data_type]]$files
    sample_numbers <- pool_data_type_map[[pool_data_type]]$sample_numbers
    lane_numbers <- pool_data_type_map[[pool_data_type]]$lane_numbers
    
    # Sort the files by both sample number and lane number
    sorted_indices <- order(sample_numbers, lane_numbers)
    sorted_file_paths <- file_paths[sorted_indices]
    
    sample_number_map <- list()
    lane_number_map <- list()
    
    # Extract unique sample and lane numbers
    for (file_path in sorted_file_paths) {
      file_name <- basename(file_path)
      sample_number <- str_extract(file_name, "(?<=AS-)\\d+")
      lane_number <- str_extract(file_name, "(?<=LR-)\\d+")
      
      # Assign unique numbers for samples and lanes starting from 1
      if (is.null(sample_number_map[[sample_number]])) {
        sample_number_map[[sample_number]] <- length(sample_number_map) + 1
      }
      if (is.null(lane_number_map[[lane_number]])) {
        lane_number_map[[lane_number]] <- length(lane_number_map) + 1
      }
    }
    
    # Store the mappings
    pool_data_type_map[[pool_data_type]]$sample_map <- sample_number_map
    pool_data_type_map[[pool_data_type]]$lane_map <- lane_number_map
  }
  
  # Create symbolic links
  for (original_path in original_file_paths) {
    pool_number_match <- regmatches(original_path, regexpr("pool\\d+", original_path))
    pool_number <- if (length(pool_number_match) > 0) pool_number_match[1] else NA
    
    data_type_match <- regmatches(original_path, regexpr("10x_multiome_\\w+_sequencing", original_path))
    data_type <- if (length(data_type_match) > 0) sub("10x_multiome_", "", sub("_sequencing", "", data_type_match[1])) else NA
    
    file_name <- basename(original_path)
    read_type <- str_extract(file_name, "_(R1|R2|I2)")
    sample_number <- str_extract(file_name, "(?<=AS-)\\d+")
    lane_number <- str_extract(file_name, "(?<=LR-)\\d+")
    
    # Log extracted information
    cat("File:", original_path, "\nExtracted info - Data type:", data_type, " Pool:", pool_number, " Sample number:", sample_number, " Lane number:", lane_number, " Read type:", read_type, "\n")
    
    if (!is.na(pool_number) && !is.na(data_type) && !is.na(read_type) && !is.na(sample_number) && !is.na(lane_number)) {
      # Use stored mappings to assign correct S and L numbers
      pool_data_type <- paste(pool_number, data_type, sep = "_")
      sample_number_id <- sprintf("S%02d", pool_data_type_map[[pool_data_type]]$sample_map[[sample_number]])
      lane_number_id <- sprintf("L%03d", pool_data_type_map[[pool_data_type]]$lane_map[[lane_number]])
      
      # Construct the new filename
      new_file_name <- paste0(pool_number, "_", sample_number_id, "_", lane_number_id, "_", substr(read_type, 2, 3), "_001.fastq.gz")
      
      # Full path for the new symbolic link
      link_dir <- file.path(link_base_dir, data_type, pool_number)
      full_new_path <- file.path(link_dir, new_file_name)
      
      # Create the directory if it doesn't exist
      dir.create(link_dir, recursive = TRUE, showWarnings = FALSE)
      
      # Creating the symbolic link
      command <- paste("ln -s", shQuote(original_path), shQuote(full_new_path))
      system(command)
      
      # Log the creation of the symlink
      cat("Created symlink:", full_new_path, "->", original_path, "\n")
    } else {
      message("Could not extract necessary information from: ", original_path)
      message("Data type:", data_type, " Pool number:", pool_number, " Read type:", read_type, " Sample number:", sample_number, " Lane number:", lane_number)
    }
  }
}

# Directory to create symbolic links
link_base_dir <- '/omics/odcf/analysis/OE0290_projects/heroes-aya_pools/AG_Thongjuea/Dataset/10x_multiomics/linked_fastqs/'

# Create symlinks for all files
create_symlink_with_rename(file_paths, link_base_dir)
