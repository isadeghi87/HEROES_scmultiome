#!/bin/bash
module load bcftools/1.9
module load htslib/1.8
 
dir=/home/i439h/projects/heroes-aya/AG_Thongjuea/Sub_project1/Result/bulkWGS/BJ-WGS/
input=/home/i439h/projects/heroes-aya-pools/temp_analysis/scRNA-Seq/Demultiplexing_results/souporcell/input_files/vcf_files

cd $dir


for i in `ls`; do
    out=${input}/$i
    cd $i/tertiary_analyses/variant_annotation/snpeff
    ### filter protein coding variants
    grep -E "(protein_coding|^#)" *snpEff.ann.vcf | grep -E '^#|rs' | grep -vE 'chr([1-9]|1[0-9]|2[0-2]|X|Y)_' > $out.protein_snp.vcf
    bcftools view -i 'QUAL>2000 & INFO/DP>100' $out.protein_snp.vcf > $out.filtered_snps.vcf
    #sed -i 's/,<NON_REF>//g' $out.high_qual.vcf
    #sed -i 's/,<NON_REF>//g' $out.filtered_snps.vcf
    ## filter only biallelic variants
    #bcftools view -m2 -M2 -v snps $out.high_qual.vcf -o $out.filtered_snps.vcf
    bgzip $out.filtered_snps.vcf
    tabix -p vcf $out.filtered_snps.vcf.gz
    rm $out.protein_snp.vcf
    cd -

done