#!/bin/bash

module load anaconda3/2021.05

module load htslib/1.8
module load bcftools/1.9
module load picard/2.25.1

chain=/home/i439h/projects/pools/AG_Thongjuea/Code/10x_multiomics/Sadeghi/Demultiplexing/references/liftover/hg19ToHg38.over.chain 
out=/home/i439h/projects/pools/AG_Thongjuea/Code/10x_multiomics/Sadeghi/Demultiplexing/liftover/vcf_files
dir=/b06x-isi/b062/g-i/HEROES-AYA/INFORM_HEROES_SNP_Genotyping
fasta=/home/i439h/projects/pools/AG_Thongjuea/Code/10x_multiomics/Sadeghi/Demultiplexing/references/hg38/genome.fa

cd $out
while IFS=$'\t' read -r col1 col2;
    
 do 
    vcf=${dir}/${col1}/snvs_${col2}_raw.vcf.gz
    vcf_out=${col1}_hg38.vcf 
    zcat $vcf | cut -f 1-10 > ./temp.vcf
   
   ## filter vcf
   picard.sh FilterVcf -I temp.vcf -O filtered.vcf
   
   ## add chr to chromosome
   awk '{if($0 !~ /^#/) print "chr"$0; else print $0}' filtered.vcf > input.vcf
   ## run liftover
   picard.sh LiftoverVcf  -I input.vcf \
        -O out.vcf \
        -C $chain \
        --REJECT ${col1}_rejected.vcf \
        -R hg38.fa \
        --ALLOW_MISSING_FIELDS_IN_HEADER true \
        --MAX_RECORDS_IN_RAM 50000 \
        --COMPRESSION_LEVEL 5
   
   rm filtered* temp.vcf input.vcf
   
   # rename header sample
   grep "CHROM" out.vcf | awk '{print $10}'> samp
   echo "$col2" > id
   paste samp id > reheader.txt

   bcftools reheader -s reheader.txt out.vcf -o $vcf_out
   bgzip $vcf_out
   tabix ${vcf_out}.gz
   rm samp id reheader.txt
done < vcf_input.txt


