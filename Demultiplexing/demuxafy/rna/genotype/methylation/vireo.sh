#!/bin/bash
set -euo pipefail

# Load necessary modules just once at the start
module load anaconda3/2021.05 singularity/3.10.5 bcftools/1.9 R/4.3.0

# Get the pool from command-line arguments
pool="${1:?Pool not specified}"
tool='vireo'

# Set project-specific directories and configurations
PROJECT_DIR="/home/i439h/projects/pools/AG_Thongjuea"
RESULTS_DIR="${PROJECT_DIR}/Result/10x_multiomics/Demultiplexing_results/demuxafy/RNA/genotype"
BAM_DIR="${PROJECT_DIR}/Result/10x_multiomics/CellRanger_ARC"
DEMUX_DIR="${PROJECT_DIR}/Software/10x_multiomics/demuxafy"
REF_DIR="${PROJECT_DIR}/Code/10x_multiomics/Sadeghi/Demultiplexing/references/hg38"
#VCF_DIR="${PROJECT_DIR}/Dataset/10x_multiomics/vcfs"
#VCF_DIR=/home/i439h/projects/pools/AG_Thongjuea/Code/10x_multiomics/Sadeghi/Demultiplexing/liftover/vcf_files
VCF_DIR=/home/i439h/projects/EMseq/results/snv_biscuit/I070-032_plasma-02-01

# Define paths
OUT_DIR="/genotype/methylation/${tool}/${pool}"
BAM="/CellRanger_ARC/$pool/outs/gex_possorted_bam.bam"
BARCODES="/CellRanger_ARC/$pool/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"
FASTA="/hg38/genome.fa"
#VCF=/vcfs/genome1K.phase3.SNP_AF5e2.hg38.vcf
#VCF=/vcf_files/I070-032_plasma-02-01/I070-032_plasma-02-01.vcf
#VCF=/vcf_files/${pool}.merged.vcf
VCF=/I070-032_plasma-02-01/I070-032_plasma-02-01.vcf.gz
THREADS=30
N=4

# Ensure output directory exists
result="${RESULTS_DIR}/methylation/$tool/$pool"
mkdir -p "$result" || { echo "Failed to create output directory: $result"; exit 1; }

# Change to the demuxafy directory
cd "$DEMUX_DIR" 

# Set up singularity bind mounts
SINGULARITY_BINDS="--bind ${RESULTS_DIR}:/genotype --bind ${BAM_DIR}:/CellRanger_ARC --bind ${REF_DIR}:/hg38  --bind ${VCF_DIR}:/I070-032_plasma-02-01"

# singularity exec $SINGULARITY_BINDS Demuxafy.sif cellsnp-lite  \
#         -b "$BARCODES" \
#         -O "$OUT_DIR" \
#         -s "$BAM"\
#         -p "$THREADS" \
#         -R $VCF \
#         --minMAF 0.1 \
#         --minCOUNT 100 \
#         --gzip \
#         --genotype

#VCF=/home/i439h/projects/EMseq/results/snv_biscuit/I070-032_plasma-02-01/I070-032_plasma-02-01.vcf.gz
# CELL="${PROJECT_DIR}/Result/10x_multiomics/Demultiplexing_results/demuxafy/RNA/genotype/methylation/${tool}/${pool}/cellSNP.cells.vcf.gz"
# out=${PROJECT_DIR}/Result/10x_multiomics/Demultiplexing_results/demuxafy/RNA/genotype/genotype/methylation/${tool}/${pool}/
# bcftools view $VCF -R $CELL -Oz -o ${out}/input.vcf.gz


singularity exec $SINGULARITY_BINDS Demuxafy.sif vireo \
    -c "$OUT_DIR" \
    -o "$OUT_DIR"\
    --callAmbientRNAs \
    -d $VCF \
    --genoTag=GT \
    -p "$THREADS" \
    --forceLearnGT \
    -N 4