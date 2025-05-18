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

pool='pool4'

# Define input files
scRNA_bam=/home/i439h/projects/pools/AG_Thongjuea/Result/10x_multiomics/CellRanger_ARC/${pool}/outs/gex_possorted_bam.bam
scATAC_bam=/home/i439h/projects/pools/AG_Thongjuea/Result/10x_multiomics/CellRanger_ARC/${pool}/outs/atac_possorted_bam.bam
genome=/home/i439h/projects/pools/AG_Thongjuea/Code/10x_multiomics/Sadeghi/Demultiplexing/references/hg38/genome.fa.fai
extract_barcodes=/home/i439h/projects/pools/AG_Thongjuea/Code/10x_multiomics/Sadeghi/merge_bams/extract_barcodes.sh

samtools merge ${pool}_merged.bam $scRNA_bam $scATAC_bam

samtools sort -o ${pool}_merged_sorted.bam ${pool}_merged.bam
samtools index ${pool}_merged_sorted.bam
