#!/bin/bash

# # Ensure the correct number of arguments are passed
# if [ "$#" -ne 1 ]; then
#     echo "Usage: $0 <PoolNumber>"
#     exit 1
# fi

# Arguments passed from the main script
PoolNumber=$1

# Load configuration
source /home/i439h/projects/heroes-aya-pools/AG_Thongjuea/Code/10x_multiomics/Sadeghi/demultiplex_pipeline/scripts/config.sh

# Validate the existence of necessary files and directories
VCF="${VCF_OUTPUT_DIR}/${PoolNumber}_filtered.vcf.gz"  # Assuming you're using the filtered VCF
if [ ! -f "$VCF" ]; then
    echo "VCF file does not exist: $VCF"
    exit 1
fi

CELL_DATA="${CELLSNP_OUTPUT_DIR}/${PoolNumber}"
if [ ! -d "$CELL_DATA" ]; then
    echo "cellSNP output directory does not exist: $CELL_DATA"
    exit 1
fi

# Define the output directory for Vireo
POOL_OUT_DIR="${VIREO_OUTPUT_DIR}/${PoolNumber}"

# Create output directory if it doesn't exist
if [ ! -d "$POOL_OUT_DIR" ]; then
        mkdir -p "$POOL_OUT_DIR"
fi

# Source necessary modules or activate Conda environment for Vireo
# Replace the following lines with the specific module load or Conda activation commands for Vireo
source ~/.bash_profile
module load bcftools/1.9
module load vcftools/0.1.16

# Run Vireo
vireo -c "$CELL_DATA" -d "$VCF" -N $n_donor --genoTag=GT --randSeed=12 -o "$POOL_OUT_DIR" \
      && echo "Vireo completed successfully for Pool $PoolNumber" \
      || { echo "Vireo failed for Pool $PoolNumber"; exit 1; }

# Check if the script execution was successful
if [ $? -eq 0 ]; then
    # Create a completion file
    touch "$OUT_DIR/vireo.completed"
    echo "vireo step completed successfully."
else
    echo "vireo step failed."
    exit 1
fi