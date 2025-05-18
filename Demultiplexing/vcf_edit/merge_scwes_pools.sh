#!/bin/bash

module load htslib/1.8
module load bcftools/1.9

cd /home/i439h/projects/pools/AG_Thongjuea/Dataset/10x_multiomics/vcfs/scWES

for i in `ls *_renamed.vcf.gz`;do tabix -p vcf $i;done
## pool 3
bcftools merge *I005_014*_renamed.vcf.gz *I010_013*_renamed.vcf.gz *I020_001*_renamed.vcf.gz *I023_037*_renamed.vcf.gz -o pool3.merged.vcf
sed -i 's/,<NON_REF>//g' pool3.merged.vcf 
bgzip pool3.merged.vcf
tabix -p vcf pool3.merged.vcf.gz

## pool 4
bcftools merge *I005_016*_renamed.vcf.gz *I010_022*_renamed.vcf.gz *I020_004*_renamed.vcf.gz *I023_047*_renamed.vcf.gz -o pool4.merged.vcf 
sed -i 's/,<NON_REF>//g' pool4.merged.vcf 
bgzip pool4.merged.vcf
tabix -p vcf pool4.merged.vcf.gz
