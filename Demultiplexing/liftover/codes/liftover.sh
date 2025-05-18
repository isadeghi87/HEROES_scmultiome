#!/bin/bash

module load anaconda3/2021.05
source activate vireo_env
module load htslib/1.8
module load bcftools/1.9
module load picard/2.25.1

chain=/home/i439h/projects/heroes-aya_pools/scRNA-Seq/Demultiplexing_results/references/liftover/hg19ToHg38.over.chain 
dir=/b06x-isi/b062/g-i/HEROES-AYA/INFORM_HEROES_SNP_Genotyping
out=/home/i439h/projects/heroes-aya/scRNA-Seq/Demultiplexing_results/liftover/vcf_files
fasta=/home/i439h/projects/heroes-aya/scRNA-Seq/CellRanger_v7.1_results/References/10x/refdata-gex-GRCh38-2020-A/fasta/genome.fa
cd $out
while IFS=$'\t' read -r col1 col2;
    
 do 
     vcf=${dir}/${col1}/snvs_${col2}_raw.vcf.gz
    vcf_out=${col1}_hg38.vcf 
    zcat $vcf | cut -f 1-8 > ./temp.vcf
   
   ## filter vcf
   picard.sh FilterVcf -I temp.vcf -O filtered.vcf
   awk '{if($0 !~ /^#/) print "chr"$0; else print $0}' filtered.vcf > input.vcf
  
    

    CrossMap.py vcf --chromid l --no-comp-allele --compress $chain input.vcf $fasta  $vcf_out
    ## compress the VCF file using bgzip and then index it using tabix
    #bcftools sort -o ${out_vcf}.gz ${out_vcf}.gz  
    #tabix -p vcf ${vcf_out}.gz
done < ${out}/vcf_input.txt


