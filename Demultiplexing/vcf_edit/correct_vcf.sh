#!/bin/bash
module load htslib/1.8
module load bcftools/1.9

cd /home/i439h/projects/pools/AG_Thongjuea/Dataset/10x_multiomics/vcfs/scWES

# Function to correct a VCF file
correct_vcf() {
    local input_vcf_gz=$1
    local corrected_vcf=${input_vcf_gz%.gz}

    # Check if the input file exists
    if [ ! -f "$input_vcf_gz" ]; then
        echo "Error: File $input_vcf_gz not found!"
        return 1
    fi

    # Remove the malformed ALT line and create a new VCF file
    zcat "$input_vcf_gz" | grep -v '##ALT=<ID=NON_REF,Description=Represents any possible alternative allele at this location>' > "$corrected_vcf"

    # Compress the corrected VCF file
    bgzip "$corrected_vcf"

    # Index the compressed VCF file
    tabix -p vcf "${corrected_vcf}.gz"

    echo "Corrected VCF file created: ${corrected_vcf}.gz"
}

# Correct pool3.merged.vcf.gz
correct_vcf "pool3.merged.vcf.gz"

# Correct pool4.merged.vcf.gz
correct_vcf "pool4.merged.vcf.gz"
