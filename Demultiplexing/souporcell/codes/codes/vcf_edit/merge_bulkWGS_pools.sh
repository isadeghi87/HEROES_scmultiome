#!/bin/bash

module load htslib/1.8
module load bcftools/1.9

dir=/home/i439h/projects/heroes-aya-pools/temp_analysis/scRNA-Seq/Demultiplexing_results/souporcell/input_files/vcf_files
cd $dir
## pool 3
bcftools merge *I005_014*filtered_snps.vcf.gz *I010_013*filtered_snps.vcf.gz *I020_001*filtered_snps.vcf.gz *I023_037*filtered_snps.vcf.gz -o pool3.merged.vcf
sed -i 's/,<NON_REF>//g' pool3.merged.vcf 
bgzip pool3.merged.vcf
tabix -p vcf pool3.merged.vcf.gz
## pool 4
bcftools merge *I005_016*filtered_snps.vcf.gz *I010_022*filtered_snps.vcf.gz *I020_004*filtered_snps.vcf.gz *I023_047*filtered_snps.vcf.gz -o pool4.merged.vcf 
sed -i 's/,<NON_REF>//g' pool4.merged.vcf 
bgzip pool4.merged.vcf
tabix -p vcf pool4.merged.vcf.gz
##### merge files with qual3000_dp100 

## pool 3
#bcftools merge *I005_014*qual3000_dp100.vcf.gz *I010_013*qual3000_dp100.vcf.gz *I020_001*qual3000_dp100.vcf.gz *I023_037*qual3000_dp100.vcf.gz -Oz -o pool3.qual3000_dp100.merged.vcf.gz 

## pool 4
#bcftools merge *I005_016*qual3000_dp100.vcf.gz *I010_022*qual3000_dp100.vcf.gz *I020_004*qual3000_dp100.vcf.gz *I023_047*qual3000_dp100.vcf.gz -Oz -o pool4.qual3000_dp100.merged.vcf.gz 

