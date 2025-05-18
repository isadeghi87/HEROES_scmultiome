#!/bin/bash
module load bcftools/1.9
module load htslib/1.8
cd /home/i439h/projects/heroes-aya-pools/temp_analysis/scRNA-Seq/Demultiplexing_results/souporcell/input_files/vcf_files

# Array of filenames to process
files=(
OE0290_HEROES-AYA_I005_014.filtered_snps.vcf.gz
OE0290_HEROES-AYA_I005_016.filtered_snps.vcf.gz
OE0290_HEROES-AYA_I010_013.filtered_snps.vcf.gz
OE0290_HEROES-AYA_I010_022.filtered_snps.vcf.gz
OE0290_HEROES-AYA_I020_001.filtered_snps.vcf.gz
OE0290_HEROES-AYA_I020_004.filtered_snps.vcf.gz
OE0290_HEROES-AYA_I023_037.filtered_snps.vcf.gz
OE0290_HEROES-AYA_I023_047.filtered_snps.vcf.gz
)

for file in "${files[@]}"
do
  output="${file/.filtered_snps.vcf.gz/.qual3000_dp100.vcf}"
  bcftools view -i 'QUAL>3000 & INFO/DP>100' $file -o $output
  bgzip $output
  tabix -p vcf "$output.gz"
done
