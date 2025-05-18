#!/bin/bash

# Enable debugging
set -x

# Load necessary modules
module load bedtools/2.24.0
module load samtools/1.5
module load parallel/20210622

# Set working directory
WD=/home/i439h/projects/pools/AG_Thongjuea/Code/10x_multiomics/Sadeghi/merge_bams/merge_dir/
cd $WD || { echo "Failed to change directory to $WD"; exit 1; }

# Define input files
scRNA_bam=/home/i439h/projects/pools/AG_Thongjuea/Result/10x_multiomics/CellRanger_ARC/pool3/outs/gex_possorted_bam.bam
scATAC_bam=/home/i439h/projects/pools/AG_Thongjuea/Result/10x_multiomics/CellRanger_ARC/pool3/outs/atac_possorted_bam.bam
genome=/home/i439h/projects/pools/AG_Thongjuea/Code/10x_multiomics/Sadeghi/Demultiplexing/references/hg38/genome.fa
extract_barcodes=/home/i439h/projects/pools/AG_Thongjuea/Code/10x_multiomics/Sadeghi/merge_bams/extract_barcodes.sh

# Check if input files exist
for file in "$scRNA_bam" "$scATAC_bam" "$genome" "$extract_barcodes"; do
  if [[ ! -f $file ]]; then
    echo "Required file not found: $file"
    exit 1
  fi
done

# Create output directory if it doesn't exist
mkdir -p $WD

# Step 1: Extract headers from BAM files
samtools view -H $scRNA_bam > scRNA_header.sam
samtools view -H $scATAC_bam > scATAC_header.sam

# Step 2: Split BAM files into smaller chunks
samtools view $scRNA_bam | split -l 100 - scRNA_chunk_
samtools view $scATAC_bam | split -l 100 - scATAC_chunk_

# Step 3: Add headers to each chunk
for chunk in scRNA_chunk_*; do
  cat scRNA_header.sam $chunk > ${chunk}_with_header.sam
  samtools view -bS ${chunk}_with_header.sam > ${chunk}.bam
  rm $chunk ${chunk}_with_header.sam
done

for chunk in scATAC_chunk_*; do
  cat scATAC_header.sam $chunk > ${chunk}_with_header.sam
  samtools view -bS ${chunk}_with_header.sam > ${chunk}.bam
  rm $chunk ${chunk}_with_header.sam
done

# Log the created chunk files
echo "Created scRNA chunk files:"
ls -lh scRNA_chunk_*.bam
echo "Created scATAC chunk files:"
ls -lh scATAC_chunk_*.bam

# Step 4: Process each chunk in parallel to append cell barcodes and convert to BED
ls scRNA_chunk_*.bam | parallel -j 4 $extract_barcodes {} {.}.bed
ls scATAC_chunk_*.bam | parallel -j 4 $extract_barcodes {} {.}.bed

# Log the generated BED files
echo "Generated scRNA BED files:"
ls -lh scRNA_chunk_*.bed
echo "Generated scATAC BED files:"
ls -lh scATAC_chunk_*.bed

# Step 5: Check if the barcodes BED files were generated
for type in scRNA scATAC; do
  if [[ ! -f ${type}_chunk_aa.bed ]]; then
    echo "Barcode BED files not generated for $type"
    exit 1
  fi
done

# Step 6: Concatenate all chunks back into single BED files
cat scRNA_chunk_*.bed > scRNA_barcodes.bed
cat scATAC_chunk_*.bed > scATAC_barcodes.bed

# Step 7: Check if concatenated BED files were created
for type in scRNA scATAC; do
  if [[ ! -f ${type}_barcodes.bed ]]; then
    echo "Concatenated BED file not found for $type"
    exit 1
  fi
done

# Step 8: Add headers to the concatenated BED files
cat scRNA_header.sam scRNA_barcodes.bed > scRNA_final.bed
cat scATAC_header.sam scATAC_barcodes.bed > scATAC_final.bed

# Step 9: Check if final BED files were created
for type in scRNA scATAC; do
  if [[ ! -f ${type}_final.bed ]]; then
    echo "Final BED file not created for $type"
    exit 1
