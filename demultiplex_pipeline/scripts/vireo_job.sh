#!/bin/bash

# Load configuration
source /home/i439h/projects/heroes-aya-pools/AG_Thongjuea/Code/10x_multiomics/Sadeghi/demultiplex_pipeline/scripts/config.sh

# Source necessary modules or activate conda environment for Vireo
source ~/.bash_profile
module load bcftools/1.9
module load vcftools/0.1.16

# Array to hold all Vireo job IDs
declare -a vireo_job_ids

# Read from the samplesheet and submit jobs for each unique pool
while read -r PoolNumber; do
    if [[ ! "$PoolNumber" =~ ^# && "$PoolNumber" != "PoolNumber" ]]; then
        BAMPath="$(awk -v pn="$PoolNumber" '$1 == pn {print $3; exit}' "$SAMPLESHEET_PATH")"
        BarcodePath="$(awk -v pn="$PoolNumber" '$1 == pn {print $5; exit}' "$SAMPLESHEET_PATH")"
        echo "Submitting Vireo for: $PoolNumber"
        job_id=$(bsub -cwd $(pwd) -q $QUEUE_NAME -R "rusage[mem=$VIREO_JOB_MEM]" -n $VIREO_JOB_CPUS -J "vireo_${PoolNumber}" \
                 "bash vireo_run.sh '$PoolNumber' '$BAMPath' '$BarcodePath'" | grep -oP '(?<=<)\d+(?=>)')
        vireo_job_ids+=($job_id)
    fi
done < <(cut -f1 "$SAMPLESHEET_PATH" | tail -n +2 | sort | uniq)

# Wait for all Vireo jobs to complete
for job_id in "${vireo_job_ids[@]}"; do
    bwait -w "done($job_id)"
done

echo "All Vireo jobs completed."
