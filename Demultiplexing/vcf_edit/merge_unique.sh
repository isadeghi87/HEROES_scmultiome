#!/bin/bash

module load htslib/1.8
module load bcftools/1.9

dir=/omics/odcf/analysis/OE0290_projects/heroes-aya_pools/AG_Thongjuea/Result/10x_multiomics/Demultiplexing_results/souporcell/input_files/vcf_files

cd $dir

# Declare an associative array with pool name as key and the associated sample identifiers as value
declare -A pools
pools=( 
    ["pool3"]="I005_014 I010_013 I020_001 I023_037"
    ["pool4"]="I005_016 I010_022 I020_004 I023_047"
    # add more pools as needed
)

# Loop through each pool and process files
for pool_name in "${!pools[@]}"; do
    # Create a list of VCF files for concatenation
    vcf_files=""
    for sample_id in ${pools[$pool_name]}; do
        vcf_files+=" ${sample_id}*unique.vcf.gz"
    done
    
    # Concatenate and filter VCFs
    cat $vcf_files > ${pool_name}.cat.vcf
    bcftools view -R ${pool_name}.cat.vcf ${pool_name}.merged.vcf.gz -o ${pool_name}.unique.vcf
    sed -i 's/,<NON_REF>//g' ${pool_name}.unique.vcf 
    cp ${pool_name}.unique.vcf ${pool_name}.unique2.vcf 
    bgzip ${pool_name}.unique.vcf
    tabix -p vcf ${pool_name}.unique.vcf.gz
done
