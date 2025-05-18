
#!/bin/bash

## load modules
module load bcftools/1.9
module load anaconda3/2021.05

## define files
input_vcf=/home/i439h/projects/pools/AG_Thongjuea/Result/10x_multiomics/Demultiplexing_results/vcf_files/pool3.merged.vcf.gz
output_dir="/home/i439h/projects/pools/AG_Thongjuea/Result/10x_multiomics/Demultiplexing_results/vcf_files/split_vcfs"
mkdir -p $output_dir

# Split VCF by chromosome
for i in {1..22} X; do
    bcftools view -r "chr${i}" $input_vcf -Oz -o $output_dir/pool3_chr${i}.vcf.gz
    bcftools index $output_dir/pool3_chr${i}.vcf.gz
done
