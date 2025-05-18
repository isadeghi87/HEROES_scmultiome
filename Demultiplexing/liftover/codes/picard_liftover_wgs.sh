#!/bin/bash
module load htslib/1.8
module load bcftools/1.9
module load picard/2.25.1


chain=/home/i439h/projects/heroes-aya-pools/temp_analysis/scRNA-Seq/Demultiplexing_results/references/liftover/hg19ToHg38.over.chain 
out=/home/i439h/projects/heroes-aya-pools/temp_analysis/scRNA-Seq/Demultiplexing_results/liftover/vcf_files/WGS
dir=/omics/odcf/project/OE0290/heroes-aya/sequencing/whole_genome_sequencing/view-by-pid/
fasta=/home/i439h/projects/heroes-aya-pools/temp_analysis/scRNA-Seq/Demultiplexing_results/references/hg38/genome.fa

cd $out
while IFS=$'\t' read -r col1 col2;
    
 do 
    vcf=$col1
    vcf_out=${col2}_hg38.vcf 
    zcat $vcf | cut -f 1-10 > ./temp.vcf
   
   ## filter vcf
   ##picard.sh FilterVcf -I temp.vcf -O filtered.vcf
   ##bcftools filter -i 'ID="^rs"' temp.vcf > rs.vcf
   ## add chr to chromosome
   awk '{if($0 !~ /^#/) print "chr"$0; else print $0}' temp.vcf > input.vcf

   ## run liftover
   picard.sh LiftoverVcf  -I input.vcf \
        -O out.vcf \
        -C $chain \
        --REJECT ${col2}_rejected.vcf \
        -R ../hg38.fa \
        --ALLOW_MISSING_FIELDS_IN_HEADER true \
        --MAX_RECORDS_IN_RAM 50000 \
        --COMPRESSION_LEVEL 5
   
   rm rs* filtered* temp.vcf input.vcf *rejected.vcf
   
   # rename header sample
   grep "CHROM" out.vcf | awk '{print $10}'> samp
   echo "$col2" > id
   paste samp id > reheader.txt

   bcftools reheader -s reheader.txt out.vcf -o $vcf_out
   bgzip $vcf_out
   tabix ${vcf_out}.gz
   rm samp id reheader.txt out.vcf*
done < vcf_input.txt


