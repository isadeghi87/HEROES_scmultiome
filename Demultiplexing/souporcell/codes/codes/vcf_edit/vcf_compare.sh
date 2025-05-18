#!/bin/bash
module load htslib/1.8
module load bcftools/1.9
vcf1=/home/i439h/projects/heroes-aya-pools/temp_analysis/scRNA-Seq/Demultiplexing_results/liftover/vcf_files/I010_013_tumor001-01_control001-01_hg38.vcf.gz
vcf2=/home/i439h/projects/heroes-aya/AG_Thongjuea/Sub_project1/Result/bulkWGS/BJ-WGS/OE0290_HEROES-AYA_I010_013/tertiary_analyses/variant_annotation/snpeff/AS-1058560-LR-69403_snpEff.ann.vcf

mkdir compare_dir
cd compare_dir
cp $vcf1 $vcf2 .
bgzip $vcf1; tabix -p vcf $vcf1.gz
bgzip $vcf2; tabix -p vcf $vcf2.gz

bcftools isec -p compare_out $vcf1 $vcf2
cd compare_out
for i in `ls *.vcf`; do
bcftools stats $i >  $i.stats.txt
done