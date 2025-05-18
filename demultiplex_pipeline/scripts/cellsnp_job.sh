#!/bin/bash

# Load configuration
source /home/i439h/projects/heroes-aya-pools/AG_Thongjuea/Code/10x_multiomics/Sadeghi/demultiplex_pipeline/scripts/config.sh

# Array to hold all cellSNP job IDs
declare -a cellsnp_job_ids

# Read from the samplesheet and submit jobs for each unique pool
while read -r PoolNumber; do
    if [[ ! "$PoolNumber" =~ ^# && "$PoolNumber" != "PoolNumber" ]]; then
        # Assuming BAMPath and BarcodePath are consistent for each pool
        BAMPath="$(awk -v pn="$PoolNumber" '$1 == pn {print $3; exit}' "$SAMPLESHEET_PATH")"
        BarcodePath="$(awk -v pn="$PoolNumber" '$1 == pn {print $5; exit}' "$SAMPLESHEET_PATH")"
        echo "Submitting cellSNP for : $PoolNumber"

        job_id=$(bsub -cwd $(pwd) -q $QUEUE_NAME -R "rusage[mem=$CELLSNP_JOB_MEM]"  -n $CELLSNP_JOB_CPUS -J "cellSNP_${PoolNumber}" \
         "bash cellsnp_run.sh '$PoolNumber' '$BAMPath' '$BarcodePath'" | grep -oP '(?<=<)\d+(?=>)')
        cellsnp_job_ids+=($job_id)
    fi
done < <(cut -f1 "$SAMPLESHEET_PATH" | tail -n +2 | sort | uniq)

# Wait for all cellSNP jobs to complete
for job_id in "${cellsnp_job_ids[@]}"; do
    bwait -w "done($job_id)"
done

echo "All cellSNP jobs completed."

# Optionally, create a completion marker file
touch "$OUT_DIR/cellsnp_all_jobs_completed.marker"
