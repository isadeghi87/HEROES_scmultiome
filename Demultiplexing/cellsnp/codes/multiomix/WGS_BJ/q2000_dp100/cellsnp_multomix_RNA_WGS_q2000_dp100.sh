#!/bin/bash
### cellsnp-lite usage for 10x multiomix 1000 genome variants

module load anaconda3/2021.05
source activate /home/i439h/projects/heroes-aya-pools/temp_analysis/tools/conda/cellsnp

# Get the pool and ID from command-line arguments
sample_number=$1
ID=$2

sample="pool${sample_number}_$ID"
pool="pool${sample_number}"
OUT_DIR="/home/i439h/projects/heroes-aya-pools/AG_Thongjuea/Result/10x_multiomics/Demultiplexing_results/cellsnp/results/multiomix/RNA/WGS_BJ/q2000_dp100/${sample}"
fasta="/home/i439h/projects/heroes-aya-pools/AG_Thongjuea/Result/10x_multiomics/Demultiplexing_results/references/hg38/genome.fa"
VCF=/home/i439h/projects/heroes-aya-pools/AG_Thongjuea/Result/10x_multiomics/Demultiplexing_results/vcf_files/$pool.merged.vcf.gz

# create results subfolder for the sample 
if [ ! -d "$OUT_DIR" ]; then
  mkdir -p "$OUT_DIR"
fi

BAM="/home/i439h/projects/heroes-aya-pools/AG_Thongjuea/Result/10x_multiomics/CellRanger_ARC/${ID}_${pool}_hg38/outs/gex_possorted_bam.bam"
BARCODE="/home/i439h/projects/heroes-aya-pools/AG_Thongjuea/Result/10x_multiomics/CellRanger_ARC/${ID}_${pool}_hg38/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"

cellsnp-lite -s $BAM -b $BARCODE -R $VCF -O $OUT_DIR -p 10 --minMAF 0.1 --minCOUNT 20 --gzip --genotype
