#!/bin/bash
# Install necessary tools
module load bcftools/1.9
module load anaconda3/2021.05

env=/home/i439h/projects/pools/AG_Thongjuea/Software/10x_multiomics/impute
results=/home/i439h/projects/pools/AG_Thongjuea/Result/10x_multiomics/imputation
mkdir -p $results

## Activate conda env
source activate $env/impute_env 

## files
ref=/home/i439h/projects/pools/AG_Thongjuea/Dataset/10x_multiomics/vcfs/genome1K.phase3.SNP_AF5e2.hg38.vcf
input=/home/i439h/projects/pools/AG_Thongjuea/Dataset/10x_multiomics/vcfs/pool3.merged.vcf


