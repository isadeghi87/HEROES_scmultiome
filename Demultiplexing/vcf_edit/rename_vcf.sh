#!/bin/bash
module load htslib/1.8
module load bcftools/1.9

set -euo pipefail

log() {
    local msg="$1"
    echo "$(date '+%Y-%m-%d %H:%M:%S') - ${msg}"
}

cd /home/i439h/projects/pools/AG_Thongjuea/Dataset/10x_multiomics/vcfs/scWES

# Function to compress, sort, index, and rename sample headers in a VCF file
rename_samples() {
    local input_vcf=$1
    local sample_id=$2
    local temp_vcf="${sample_id}_temp.vcf.gz"
    local sorted_vcf="${sample_id}_sorted.vcf.gz"
    local renamed_vcf="${sample_id}_renamed.vcf.gz"

    log "Processing ${input_vcf} for sample ${sample_id}"

    # Check if input VCF file exists
    if [[ ! -f "${input_vcf}" ]]; then
        log "Error: File ${input_vcf} does not exist"
        return 1
    fi

    # Compress the VCF file with bgzip
    log "Compressing ${input_vcf}..."
    if ! bgzip -c "${input_vcf}" > "${temp_vcf}"; then
        log "Error: Failed to compress ${input_vcf}"
        return 1
    fi

    # Sort the compressed VCF file
    log "Sorting ${temp_vcf}..."
    if ! bcftools sort "${temp_vcf}" -Oz -o "${sorted_vcf}"; then
        log "Error: Failed to sort ${temp_vcf}"
        return 1
    fi

    # Index the sorted VCF file
    log "Indexing ${sorted_vcf}..."
    if ! bcftools index "${sorted_vcf}"; then
        log "Error: Failed to index ${sorted_vcf}"
        return 1
    fi

    # Rename the sample in the VCF
    log "Renaming sample in ${sorted_vcf}..."
    if ! bcftools reheader -s <(echo "${sample_id}") -o "${renamed_vcf}" "${sorted_vcf}"; then
        log "Error: Failed to rename ${sorted_vcf}"
        return 1
    fi

    # Index the renamed VCF file
    log "Indexing ${renamed_vcf}..."
    if ! bcftools index "${renamed_vcf}"; then
        log "Error: Failed to index ${renamed_vcf}"
        return 1
    fi

    log "Successfully processed ${input_vcf} for sample ${sample_id}"
}


# To process multiple files, you can loop over them like this:
vcf_files=("I005_014_merged_scWES.vcf" "I010_013_merged_scWES.vcf" "I020_001_merged_scWES.vcf" "I023_037_merged_scWES.vcf"
"I005_016_merged_scWES.vcf" "I010_022_merged_scWES.vcf"  "I020_004_merged_scWES.vcf"  "I023_047_merged_scWES.vcf")

for vcf in "${vcf_files[@]}"; do
    sample_id=$(basename "${vcf}" | cut -d'_' -f1-2)
    rename_samples "${vcf}" "${sample_id}"
done
