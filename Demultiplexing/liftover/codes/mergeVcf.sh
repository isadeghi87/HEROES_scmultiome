#!/bin/bash

## merge multiple VCF files into one file
module load bcftools/1.9
module load vcftools/0.1.16
module load htslib/1.8

# set filters
#MAF=0.1
#MISS=0.1
#QUAL=30
#MIN_DEPTH=10
#MAX_DEPTH=50

## PIDs for each pool
dir=/home/i439h/projects/heroes-aya/scRNA-Seq/Demultiplexing_results/liftover/vcf_files 
cd $dir


for pool in Pool1 Pool2 Pool3 Pool4;
do
    echo $pool
    input=$(grep -e "$pool" sampleID.tsv | cut  -f3 | paste -s -d" ")
    output=$pool.vcf.gz
    
    echo $input

    #vcf-merge --intersect $input | bgzip -c > $output
    bcftools merge $input -Oz -o out.vcf.gz


    ## keep only main chromosoes
    zgrep -vE 'chr([1-9]|1[0-9]|2[0-2]|X|Y)_' out.vcf.gz | gzip > temp.vcf.gz 
    
    # perform the filtering variants
    DP=10
    QUAL=50

    #bcftools view -v snps temp.vcf.gz -Oz -o snp.vcf.gz
    #bcftools view -i "QUAL>= $QUAL || INFO/DP > $DP " temp.vcf.gz -Oz -o $output
    bcftools filter -i "QUAL>= $QUAL" temp.vcf.gz -Oz -o $output
    tabix -p vcf $output
    rm out.vcf.gz temp.vcf.gz 
done



