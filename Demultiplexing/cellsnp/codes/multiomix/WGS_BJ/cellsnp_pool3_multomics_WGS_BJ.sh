#!/bin/bash
### cellsnp-lite usage
module load anaconda3/2021.05
source activate /home/i439h/projects/heroes-aya-pools/temp_analysis/tools/conda/cellsnp
sample='pool3_1016715'
ID='1016715'
pool='pool3'
OUT_DIR=/home/i439h/projects/heroes-aya-pools/temp_analysis/scRNA-Seq/Demultiplexing_results/cellsnp/results/${sample}/WGS_BJ
fasta=/home/i439h/projects/heroes-aya-pools/temp_analysis/scRNA-Seq/Demultiplexing_results/references/hg38/genome.fa

# create results subfolder for the sample 
if [ ! -d "$OUT_DIR" ]; then
  mkdir -p "$OUT_DIR"
fi

BAM=/home/i439h/projects/heroes-aya-pools/temp_analysis/scMulti-omics/CellRanger_arc/${ID}_${pool}_hg38/outs/gex_possorted_bam.bam
BARCODE=/home/i439h/projects/heroes-aya-pools/temp_analysis/scMulti-omics/CellRanger_arc/${ID}_${pool}_hg38/outs/filtered_feature_bc_matrix/barcodes.tsv.gz
VCF=/home/i439h/projects/heroes-aya-pools/temp_analysis/scRNA-Seq/Demultiplexing_results/souporcell/input_files/vcf_files/$pool.qual.vcf

cellsnp-lite -s $BAM -b $BARCODE -R $VCF -O $OUT_DIR --UMItag UB -p 22 --minMAF 0.1 --minCOUNT 20 --gzip --genotype

