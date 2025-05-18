#!/bin/bash

# Load the necessary modules
module load htslib/1.8
module load bcftools/1.9
cd /omics/odcf/analysis/OE0290_projects/heroes-aya_pools/AG_Thongjuea/Result/10x_multiomics/Demultiplexing_results/souporcell/input_files/vcf_files
# Define an associative array where keys are sample identifiers and values are pool files
declare -A samples
samples=(
    ["I005_014"]="pool3.merged2.vcf.gz"
    ["I010_013"]="pool3.merged2.vcf.gz"
    ["I020_001"]="pool3.merged2.vcf.gz"
    ["I023_037"]="pool3.merged2.vcf.gz"
    ["I005_016"]="pool4.merged2.vcf.gz"
    ["I010_022"]="pool4.merged2.vcf.gz"
    ["I020_004"]="pool4.merged2.vcf.gz"
    ["I023_047"]="pool4.merged2.vcf.gz"
)

# Loop through each sample and run bcftools isec
for sample in "${!samples[@]}"; do
    output_vcf="${sample}_unique.vcf.gz"
    pool_vcf="${samples[$sample]}"
    input_vcf="OE0290_HEROES-AYA_${sample}.filtered_snps.vcf.gz"

    # Check if the input and pool VCF files exist before running bcftools isec
    if [[ -f ${input_vcf} ]] && [[ -f ${pool_vcf} ]]; then
        bcftools isec -n=1 -Oz -o ${output_vcf} ${pool_vcf} ${input_vcf}
    else
        echo "Input VCF file ${input_vcf} or pool VCF file ${pool_vcf} not found, skipping ${sample}"
    fi
done
