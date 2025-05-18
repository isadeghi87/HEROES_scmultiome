#!/bin/bash
set -euo pipefail

# Load necessary modules just once at the start
module load anaconda3/2021.05 singularity/3.10.5 bcftools/1.9 R/4.3.0

# Get the pool from command-line arguments
pool="${1:?Pool not specified}"
tool='souporcell'

# Set project-specific directories and configurations
PROJECT_DIR="/home/i439h/projects/pools/AG_Thongjuea"
RESULTS_DIR="${PROJECT_DIR}/Result/10x_multiomics/Demultiplexing_results/demuxafy/RNA/genotype/scWES"
BAM_DIR="${PROJECT_DIR}/Result/10x_multiomics/CellRanger_ARC"
DEMUX_DIR="${PROJECT_DIR}/Software/10x_multiomics/demuxafy"
REF_DIR="${PROJECT_DIR}/Code/10x_multiomics/Sadeghi/Demultiplexing/references/hg38"
VCF_DIR="${PROJECT_DIR}/Dataset/10x_multiomics/vcfs/scWES/"

# Define paths
OUT_DIR="/genotype/scWES/${tool}/${pool}"
BAM="/CellRanger_ARC/$pool/outs/gex_possorted_bam.bam"
BARCODES="/CellRanger_ARC/$pool/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"
FASTA="/hg38/genome.fa"
VCF=/scWES/${pool}.merged.vcf
THREADS=20
N=4

# Ensure output directory exists
result="${RESULTS_DIR}/$tool/$pool"
mkdir -p "$result" 

# Change to the demuxafy directory
cd "$DEMUX_DIR" 

# Set up singularity bind mounts
SINGULARITY_BINDS="--bind ${RESULTS_DIR}:/genotype/scWES --bind ${BAM_DIR}:/CellRanger_ARC --bind ${REF_DIR}:/hg38 --bind ${VCF_DIR}:/scWES"

# Execute the souporcell pipeline using singularity with appropriate bindings
singularity exec $SINGULARITY_BINDS Demuxafy.sif souporcell_pipeline.py \
                -i "$BAM" \
                -b "$BARCODES" \
                -f "$FASTA" \
                -t "$THREADS" \
                -o "$OUT_DIR" \
                -k "$N" \
               --known_genotypes $VCF

#Souporcell Summary
singularity exec $SINGULARITY_BINDS Demuxafy.sif bash souporcell_summary.sh $OUT_DIR/clusters.tsv > /omics/odcf/analysis/OE0290_projects/heroes-aya_pools/AG_Thongjuea/Result/10x_multiomics/Demultiplexing_results/demuxafy/RNA/genotype/scWES/$tool/$pool/${tool}_summary.tsv

echo "Pool: $pool"

Rscript souporcell_assign.R "$pool"