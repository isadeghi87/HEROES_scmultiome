#!/bin/bash
## This script extracts reads from a pooled BAM file based on cell barcodes associated with a specific sample
## and then sorts and indexes the resulting BAM files.

# Load necessary modules
module load samtools/1.5

# Define paths and filenames
WORK_DIR="/home/i439h/projects/pools/AG_Thongjuea/Result/10x_multiomics/CellRanger_ARC/pool12/outs"
BARCODES_FILE="/home/i439h/projects/pools/AG_Thongjuea/Result/10x_multiomics/integrative_analysis/tables/I070_032_barcodes.tsv"
RNA_BAM="gex_possorted_bam.bam"
ATAC_BAM="atac_possorted_bam.bam"
SAMPLE="I070_032"
OUT_ATAC="${SAMPLE}_atac.bam"
OUT_RNA="${SAMPLE}_rna.bam"
SORTED_ATAC="${SAMPLE}_atac_sorted.bam"
SORTED_RNA="${SAMPLE}_rna_sorted.bam"

# Change to the working directory
cd $WORK_DIR

# Debug: Check if barcode file is non-empty and has the correct format
if [ ! -s $BARCODES_FILE ]; then
    echo "Error: Barcode file is empty or not found."
    exit 1
fi

# Extract RNA reads based on cell barcodes
samtools view -h $RNA_BAM | \
awk 'NR==FNR {barcodes["CB:Z:" $1]=1; next} /^@/ || ($0 ~ /CB:Z:/ && $0 ~ /CB:Z:[A-Za-z0-9-]+/ && substr($0, match($0, /CB:Z:[A-Za-z0-9-]+/), RLENGTH) in barcodes)' $BARCODES_FILE - | \
samtools view -b -o $OUT_RNA -

# Debug: Count number of reads in output
echo "RNA reads extracted: $(samtools view -c $OUT_RNA)"

# Sort and index the RNA BAM file
samtools sort -o $SORTED_RNA $OUT_RNA
samtools index $SORTED_RNA

# Extract ATAC reads based on cell barcodes
samtools view -h $ATAC_BAM | \
awk 'NR==FNR {barcodes["CB:Z:" $1]=1; next} /^@/ || ($0 ~ /CB:Z:/ && $0 ~ /CB:Z:[A-Za-z0-9-]+/ && substr($0, match($0, /CB:Z:[A-Za-z0-9-]+/), RLENGTH) in barcodes)' $BARCODES_FILE - | \
samtools view -b -o $OUT_ATAC -

# Debug: Count number of reads in output
echo "ATAC reads extracted: $(samtools view -c $OUT_ATAC)"

# Sort and index the ATAC BAM file
samtools sort -o $SORTED_ATAC $OUT_ATAC
samtools index $SORTED_ATAC
