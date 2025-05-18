#!/bin/bash

vcf=/home/i439h/projects/heroes-aya-pools/temp_analysis/scRNA-Seq/Demultiplexing_results/souporcell/input_files/vcf_files/pool3.merged2.vcf.gz
bam=/home/i439h/projects/heroes-aya-pools/AG_Thongjuea/Result/10x_multiomics/CellRanger_ARC/1016715_pool3_hg38/outs/gex_possorted_bam.bam
barcode=/home/i439h/projects/heroes-aya-pools/AG_Thongjuea/Result/10x_multiomics/CellRanger_ARC/1016715_pool3_hg38/outs/filtered_feature_bc_matrix/barcodes.tsv.gz
out=/home/i439h/projects/heroes-aya-pools/AG_Thongjuea/Code/10x_multiomics/Sadeghi/demultiplex_pipeline/test
# Load neccassary modules
module load anaconda3/2021.05
source activate /home/i439h/projects/heroes-aya-pools/temp_analysis/tools/conda/cellsnp

cellsnp-lite -s $bam -b $barcode -R $vcf -O $out \
             -p 10 --minMAF 0.1 --minCOUNT 20 --gzip --genotype