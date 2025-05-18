#!/bin/bash
set -euo pipefail

# Load necessary modules just once at the start
module load anaconda3/2021.05 singularity/3.10.5 bcftools/1.9 R/4.3.0

# Get the pool from command-line arguments
pool="${1:?Pool not specified}"
tool='souporcell'

# Set project-specific directories and configurations
PROJECT_DIR="/home/i439h/projects/pools/AG_Thongjuea"
RESULTS_DIR="${PROJECT_DIR}/Result/10x_multiomics/Demultiplexing_results/demuxafy/merged/genotype"
BAM_DIR="${PROJECT_DIR}/Code/10x_multiomics/Sadeghi/merge_bams/merge_dir/"
BAR_DIR="${PROJECT_DIR}/Result/10x_multiomics/CellRanger_ARC"
DEMUX_DIR="${PROJECT_DIR}/Software/10x_multiomics/demuxafy"
REF_DIR="${PROJECT_DIR}/Code/10x_multiomics/Sadeghi/Demultiplexing/references/hg38"
VCF_DIR="${PROJECT_DIR}/Result/10x_multiomics/Demultiplexing_results/vcf_files"

# Define paths
OUT_DIR="/genotype/${tool}/${pool}"
BAM=/merge_dir/${pool}_merged_sorted.bam
BARCODES="/CellRanger_ARC/$pool/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"
FASTA="/hg38/genome.fa"
VCF=/vcf_files/${pool}.merged.vcf
THREADS=20
N=4

# Ensure output directory exists
result="${RESULTS_DIR}/$tool/$pool"
mkdir -p "$result" 

# Change to the demuxafy directory
cd "$DEMUX_DIR" 

# Set up singularity bind mounts
SINGULARITY_BINDS="--bind ${RESULTS_DIR}:/genotype --bind ${BAR_DIR}:/CellRanger_ARC  --bind ${BAM_DIR}:/merge_dir --bind ${REF_DIR}:/hg38 --bind ${VCF_DIR}:/vcf_files"

# Execute the souporcell pipeline using singularity with appropriate bindings
singularity exec $SINGULARITY_BINDS Demuxafy.sif souporcell_pipeline.py \
                -i "$BAM" \
                -b "$BARCODES" \
                -f "$FASTA" \
                -t "$THREADS" \
                -o "$OUT_DIR" \
                -k "$N" \
                 --known_genotypes $VCF \
                --no_umi True \
                --skip_remap True

#Souporcell Summary
singularity exec $SINGULARITY_BINDS Demuxafy.sif bash souporcell_summary.sh $OUT_DIR/clusters.tsv > /omics/odcf/analysis/OE0290_projects/heroes-aya_pools/AG_Thongjuea/Result/10x_multiomics/Demultiplexing_results/demuxafy/merged/genotype/$tool/$pool/${tool}_summary.tsv

#Correlating Cluster to Donor Reference SNP Genotypes
singularity exec --bind $PROJECT_DIR/Result/10x_multiomics/Demultiplexing_results/vcf_files:/vcf_files \
        --bind $PROJECT_DIR/Result/10x_multiomics/Demultiplexing_results/demuxafy/merged/genotype/$tool:/$tool Demuxafy.sif Assign_Indiv_by_Geno.R \
        -r $VCF -c $OUT_DIR/cluster_genotypes.vcf -o $OUT_DIR

 echo "Pool: $pool"

# #Rscript souporcell_assign.R "$pool"
