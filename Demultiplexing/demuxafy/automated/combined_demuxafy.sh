#!/bin/bash
set -euo pipefail

# Load necessary modules just once at the start
module load anaconda3/2021.05 singularity/3.10.5 bcftools/1.9 R/4.3.0

# Get the pool from command-line arguments
pool="${1:?Pool not specified}"

# Set project-specific directories and configurations
PROJECT_DIR="/home/i439h/projects/pools/AG_Thongjuea"
RESULTS_DIR_BASE="${PROJECT_DIR}/Result/10x_multiomics/Demultiplexing_results/demuxafy"
BAM_DIR="${PROJECT_DIR}/Result/10x_multiomics/CellRanger_ARC"
DEMUX_DIR="${PROJECT_DIR}/Software/10x_multiomics/demuxafy"
REF_DIR="${PROJECT_DIR}/Code/10x_multiomics/Sadeghi/Demultiplexing/references/hg38"
VCF_DIR="${PROJECT_DIR}/Dataset/10x_multiomics/vcfs"
FASTA="/hg38/genome.fa"
VCF="/vcfs/genome1K.phase3.SNP_AF5e2.hg38.vcf"

# Set up singularity bind mounts
SINGULARITY_BINDS="--bind ${RESULTS_DIR_BASE}:/no_genotype --bind ${BAM_DIR}:/CellRanger_ARC --bind ${REF_DIR}:/hg38 --bind ${VCF_DIR}:/vcfs"

run_souporcell() {
  local modality=$1
  local threads=20
  local n=4
  local tool='souporcell'
  local results_dir="${RESULTS_DIR_BASE}/${modality}/no_genotype"
  local bam
  local barcodes

  if [ "$modality" == "RNA" ]; then
    bam="/CellRanger_ARC/$pool/outs/gex_possorted_bam.bam"
    barcodes="/CellRanger_ARC/$pool/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"
  elif [ "$modality" == "ATAC" ]; then
    bam="/CellRanger_ARC/$pool/outs/atac_possorted_bam.bam"
    barcodes="/CellRanger_ARC/$pool/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"
  fi

  local out_dir="/no_genotype/${tool}/${pool}"
  local result="${results_dir}/$tool/$pool"
  mkdir -p "$result"

  cd "$DEMUX_DIR"

  singularity exec $SINGULARITY_BINDS Demuxafy.sif souporcell_pipeline.py \
                  -i "$bam" \
                  -b "$barcodes" \
                  -f "$FASTA" \
                  -t "$threads" \
                  -o "$out_dir" \
                  -k "$n" \
                  --common_variants $VCF \
                  --no_umi True \
                  --skip_remap True

  singularity exec $SINGULARITY_BINDS Demuxafy.sif bash souporcell_summary.sh $out_dir/clusters.tsv > $results_dir/$tool/$pool/${tool}_summary.tsv
}

run_vireo() {
  local modality=$1
  local threads=30
  local tool='vireo'
  local results_dir="${RESULTS_DIR_BASE}/${modality}/no_genotype"
  local bam
  local barcodes

  if [ "$modality" == "RNA" ]; then
    bam="/CellRanger_ARC/$pool/outs/gex_possorted_bam.bam"
    barcodes="/CellRanger_ARC/$pool/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"
  elif [ "$modality" == "ATAC" ]; then
    bam="/CellRanger_ARC/$pool/outs/atac_possorted_bam.bam"
    barcodes="/CellRanger_ARC/$pool/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"
  fi

  local out_dir="/no_genotype/${tool}/${pool}"
  local result="${results_dir}/$tool/$pool"
  mkdir -p "$result"

  cd "$DEMUX_DIR"

  singularity exec $SINGULARITY_BINDS Demuxafy.sif cellsnp-lite \
          -b "$barcodes" \
          -O "$out_dir" \
          -s "$bam" \
          -p "$threads" \
          -R $VCF \
          --minMAF 0.1 \
          --minCOUNT 20 \
          --gzip \
          --genotype \
          --cellTAG CB \
          --UMItag None

  singularity exec $SINGULARITY_BINDS Demuxafy.sif vireo \
      -c "$out_dir" \
      -o "$out_dir" \
      --callAmbientRNAs \
      --genoTag=GT \
      -p "$threads" \
      -N 4

  # Generate Vireo Summary (assuming there is a similar script for Vireo summary)
  singularity exec $SINGULARITY_BINDS Demuxafy.sif bash vireo_summary.sh $out_dir > $results_dir/$tool/$pool/${tool}_summary.tsv
}

combine_results() {
  local modality=$1
  local tool='combined'
  local results_dir="${RESULTS_DIR_BASE}/${modality}/no_genotype/$pool"
  local souporcell_outdir="${results_dir}/souporcell/$pool"
  local vireo_outdir="${results_dir}/vireo/$pool"
  local out_dir="${results_dir}/combined"

  mkdir -p "$out_dir"

  singularity exec $SINGULARITY_BINDS Demuxafy.sif Combine_Results.R \
    -o $out_dir/combined_results.tsv \
    --souporcell $souporcell_outdir \
    --vireo $vireo_outdir \
    --method "MajoritySinglet"
}

# Run all analyses in parallel
run_souporcell "RNA" &
run_souporcell "ATAC" &
run_vireo "RNA" &
run_vireo "ATAC" &

# Wait for all background jobs to finish
wait

# Combine results for each modality
combine_results "RNA"
combine_results "ATAC"

echo "All analyses complete for pool $pool."
