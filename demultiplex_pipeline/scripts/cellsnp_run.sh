#!/bin/bash

# Load configuration
source  /home/i439h/projects/heroes-aya-pools/AG_Thongjuea/Code/10x_multiomics/Sadeghi/demultiplex_pipeline/scripts/config.sh

# Load neccassary modules
module load anaconda3/2021.05
source activate /home/i439h/projects/heroes-aya-pools/temp_analysis/tools/conda/cellsnp
# Validate parameters
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <PoolNumber> <BAMPath> <BarcodePath>"
    exit 1
fi

PoolNumber=$1
BAMPath=$2
BarcodePath=$3

# Verify input files exist
if [ ! -f "$BAMPath" ]; then
    echo "BAM file does not exist: $BAMPath"
    exit 1
fi

if [ ! -f "$BarcodePath" ]; then
    echo "Barcode file does not exist: $BarcodePath"
    exit 1
fi

# Define the input VCF file (merged VCF for the pool)
VCF="${VCF_OUTPUT_DIR}/${PoolNumber}_filtered.vcf.gz"
if [ ! -f "$VCF" ]; then
    echo "VCF file does not exist: $VCF"
    exit 1
fi

echo "Running cellSNP for pool: $PoolNumber"
echo "BAMPath: $BAMPath"
echo "BarcodePath: $BarcodePath"
echo "VCF: $VCF"

# Define the output directory for cellSNP
POOL_OUT_DIR="${CELLSNP_OUTPUT_DIR}/${PoolNumber}"

# Create output directory if it doesn't exist
if [ ! -d "$POOL_OUT_DIR" ]; then
        mkdir -p "$POOL_OUT_DIR"
fi
# Run cellSNP
echo "Running cellSNP for pool: $PoolNumber"
if cellsnp-lite -s "$BAMPath" -b "$BarcodePath" -R "$VCF" -O "$POOL_OUT_DIR" \
   -p 10 --minMAF 0.1 --minCOUNT 20 --gzip --genotype; then
    echo "cellSNP-lite completed successfully for pool: $PoolNumber."
else
    echo "cellSNP-lite failed for pool: $PoolNumber."
    exit 1
fi


# Verify the output
output_vcf="$POOL_OUT_DIR/cellsnp.base.vcf.gz"
if [ -f "$output_vcf" ]; then
    # Check for the presence of the VCF header
    if ! zcat "$output_vcf" | grep -qm1 "^#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"; then
        echo "VCF header is missing in the output file: $output_vcf"
        exit 1
    fi

    # Check for the presence of at least one variant entry
    if ! zcat "$output_vcf" | grep -qm1 "^chr"; then
        echo "No variant data found in the output file: $output_vcf"
        exit 1
    fi

    echo "cellSNP output validation passed for pool: $PoolNumber."
else
    echo "cellSNP output file does not exist: $output_vcf"
    exit 1
fi

# Check if the script execution was successful
if [ $? -eq 0 ]; then
    # Create a completion file
    touch "$OUT_DIR/cellsnp.completed"
    echo "cellsnp step completed successfully."
else
    echo "cellsnp step failed."
    exit 1
fi