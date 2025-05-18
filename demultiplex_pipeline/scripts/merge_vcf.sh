#!/bin/bash

# Load necessary modules
module load bcftools/1.9
module load htslib/1.8

# Load configuration
source  /home/i439h/projects/heroes-aya-pools/AG_Thongjuea/Code/10x_multiomics/Sadeghi/demultiplex_pipeline/scripts/config.sh

merge_and_filter_vcf() {
     local pool_id="$1"
    shift # Remove the first argument, which is pool_id
    local vcf_files=("$@") # Capture the remaining arguments as VCF files
    
    local merged_vcf="${VCF_OUTPUT_DIR}/merged_${pool_id}.vcf"
    local filtered_vcf="${VCF_OUTPUT_DIR}/${pool_id}_filtered.vcf"
   
   # Check if filtered VCF already exists to skip processing
    if [[ -f "$filtered_vcf.gz" ]]; then
        echo "$filtered_vcf.gz already exists, skipping..."
        return 0
    fi
   
    echo "Merging files for pool: $pool_id"
    printf '%s\n' "Merging the following VCF files for pool: $pool_id" "${vcf_files[@]}"

    bcftools merge "${vcf_files[@]}" -o "$merged_vcf"
    [ $? -eq 0 ] && echo "Merging completed: $merged_vcf" || { echo "Merging failed for $pool_id"; return 1; }
    

    # Temporary file for storing rs filtered VCF
    local rs_vcf=$(mktemp "${VCF_OUTPUT_DIR}/rs_${pool_id}_XXXXXX.vcf")

    # Filter VCF for rs IDs and remove non-biallelic variants, then filter by quality and depth
    echo "Keep variants with rs ID and from main chromosomes"
    bcftools view "$merged_vcf" | \
    grep -E '^#|rs' | \
    grep -vE 'chr([1-9]|1[0-9]|2[0-2]|X|Y)_' > "$rs_vcf"
    sed -i 's/,<NON_REF>//g' "$rs_vcf"

    # Perform filtering for biallelic variants and quality-depth criteria
    echo "Filtering biallelic snps"
    bcftools view -m2 -M2 -v snps "$rs_vcf" -o "${rs_vcf%.vcf}_biallelic.vcf"

    echo "Filtering low quality snps"
    bcftools view -i 'QUAL>100 & INFO/DP>20' "${rs_vcf%.vcf}_biallelic.vcf" -o "$filtered_vcf"

    echo "Compressing and indexing filtered VCF for pool: $pool_id"
    bgzip -f "$filtered_vcf"
    tabix -p vcf "$filtered_vcf.gz"

    # Final check to ensure the filtered file exists and is not empty
    if [ -s "$filtered_vcf.gz" ]; then
        echo "Final output ready: $filtered_vcf.gz"
    else
        echo "Final output file missing or empty for $pool_id"
        return 1
    fi

    # Cleanup temporary files
    echo "Cleanup temporary files"
    rm "$rs_vcf" "${rs_vcf%.vcf}_biallelic.vcf"
}



declare -A pool_vcf_files

# Read from the samplesheet and prepare VCF files for merging
while IFS=$'\t' read -r PoolNumber Sample BAMPath VCFPath BarcodePath; do
    if [[ "$PoolNumber" != "PoolNumber" ]]; then
         reheadered_vcf="${VCF_OUTPUT_DIR}/${Sample}.reheadered.vcf.gz"
        # Ensure this matches the output naming convention from your reheadering step
        if [[ -f "$reheadered_vcf" ]]; then
            pool_vcf_files["$PoolNumber"]+="$reheadered_vcf "
        else
            echo "Warning: Expected VCF file not found - $reheadered_vcf"
        fi
    fi
done < "$SAMPLESHEET_PATH"

# Process VCF files for each pool
for pool in "${!pool_vcf_files[@]}"; do
    # Convert string to array to handle file paths correctly
    IFS=' ' read -r -a vcf_files_for_pool <<< "${pool_vcf_files[$pool]}"
    if [ ${#vcf_files_for_pool[@]} -gt 0 ]; then
        merge_and_filter_vcf "$pool" "${vcf_files_for_pool[@]}"
    else
        echo "No VCF files to process for pool: $pool"
    fi
done

