#!/bin/bash

# Load configuration
module load bcftools/1.9
module load htslib/1.8

source  /home/i439h/projects/heroes-aya-pools/AG_Thongjuea/Code/10x_multiomics/Sadeghi/demultiplex_pipeline/scripts/config.sh


# Create directories if they don't exist
for dir in "$VCF_OUTPUT_DIR" "$CELLSNP_OUTPUT_DIR" "$VIREO_OUTPUT_DIR" "$SOUPORCELL_OUTPUT_DIR" "$SAMPLE_LIST_DIR"; do
    if [ ! -d "$dir" ]; then
        mkdir -p "$dir"
    fi
done

# Define maximum number of parallel jobs
max_jobs=4
current_jobs=0

# Function to process a single VCF file
process_vcf() {
    PoolNumber="$1"
    Sample="$2"
    BAMPath="$3"
    VCFPath="$4"
    BarcodePath="$5"

    reheadered_vcf="${VCF_OUTPUT_DIR}/${Sample}.reheadered.vcf"

    if [ -f "${reheadered_vcf}.gz" ]; then
        echo "File ${reheadered_vcf}.gz already exists, skipping..."
        return 0
    fi

    echo "$Sample: Reheadering and preparing VCF"
    bcftools reheader -s <(echo "$Sample") -o "$reheadered_vcf" "$VCFPath"
    bgzip -f "$reheadered_vcf"
    tabix -p vcf -f "${reheadered_vcf}.gz"
}

# Read from the samplesheet and process VCF files in parallel
while IFS=$'\t' read -r PoolNumber Sample BAMPath VCFPath BarcodePath; do
    if [ "$PoolNumber" != "PoolNumber" ]; then
        process_vcf "$PoolNumber" "$Sample" "$BAMPath" "$VCFPath" "$BarcodePath" &

        ((current_jobs++))
        if [ "$current_jobs" -ge "$max_jobs" ]; then
            wait # Wait for all background jobs to finish
            current_jobs=0
        fi
    fi
done < "$SAMPLESHEET_PATH"
wait # Wait for the last batch of background jobs to finish

echo "Filtering VCF process completed for all samples."
touch "$OUT_DIR/filter_vcf.completed"