#!/bin/bash

module load htslib/1.8
module load bcftools/1.9

dir=/home/i439h/projects/pools/AG_Thongjuea/Code/10x_multiomics/Sadeghi/Demultiplexing/liftover/vcf_files
cd $dir

# ## pool 3
# bcftools merge *I005_014*filtered_snps.vcf.gz *I010_013*filtered_snps.vcf.gz *I020_001*filtered_snps.vcf.gz *I023_037*filtered_snps.vcf.gz -o pool3.merged.vcf
# sed -i 's/,<NON_REF>//g' pool3.merged.vcf 
# bgzip pool3.merged.vcf
# tabix -p vcf pool3.merged.vcf.gz
# cp pool3.merged.vcf.gz pool3.merged2.vcf.gz
# gunzip pool3.merged2.vcf.gz

# #reheader
# bcftools query -l pool3.merged2.vcf > pool3_old.txt
# paste pool3_id.txt pool3_old.txt >reheader_pool3.txt 
# bcftools reheader -s reheader_pool3.txt -o pool3.merged.vcf pool3.merged2.vcf

# ## pool 4
# bcftools merge *I005_016*filtered_snps.vcf.gz *I010_022*filtered_snps.vcf.gz *I020_004*filtered_snps.vcf.gz *I023_047*filtered_snps.vcf.gz -o pool4.merged.vcf 
# sed -i 's/,<NON_REF>//g' pool4.merged.vcf 
# bgzip pool4.merged.vcf
# tabix -p vcf pool4.merged.vcf.gz
# cp pool4.merged.vcf.gz pool4.merged2.vcf.gz
# gunzip pool4.merged2.vcf.gz

# #reheader
# bcftools query -l pool4.merged2.vcf > pool4_old.txt
# paste pool4_id.txt pool4_old.txt >reheader_pool4.txt 
# bcftools reheader -s reheader_pool3.txt -o pool4.merged.vcf pool4.merged2.vcf


## pool 14
# bcftools merge I074_025*hg38.vcf.gz I007_072*hg38.vcf.gz -o pool14.merged.vcf 
# sed -i 's/,<NON_REF>//g' pool14.merged.vcf 
# bgzip pool14.merged.vcf
# tabix -p vcf pool14.merged.vcf.gz
# cp pool14.merged.vcf.gz pool14.merged2.vcf.gz
# gunzip pool14.merged2.vcf.gz

# #reheader
# bcftools query -l pool4.merged2.vcf > pool4_old.txt
# paste pool4_id.txt pool4_old.txt >reheader_pool4.txt 
# bcftools reheader -s reheader_pool3.txt -o pool14.merged.vcf pool4.merged2.vcf


## pool 13
bcftools merge I007_057*hg38.vcf.gz I052_006*hg38.vcf.gz -o pool13.merged.vcf 
sed -i 's/,<NON_REF>//g' pool13.merged.vcf 
bgzip pool13.merged.vcf
tabix -p vcf pool13.merged.vcf.gz
cp pool13.merged.vcf.gz pool13.merged2.vcf.gz
gunzip pool13.merged2.vcf.gz