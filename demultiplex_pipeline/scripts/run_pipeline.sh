#!/bin/bash

# Load configuration
source  /home/i439h/projects/heroes-aya-pools/AG_Thongjuea/Code/10x_multiomics/Sadeghi/demultiplex_pipeline/scripts/config.sh

# Logging setup
LOG_DIR="$OUT_DIR/logs"
mkdir -p "$LOG_DIR"
LOG_FILE="$LOG_DIR/pipeline_run_$(date +%Y%m%d_%H%M%S).log"

echo "Starting the pipeline execution..." | tee -a "$LOG_FILE"


# Function to submit a job and wait for its completion
submit_and_wait() {
    local script=$1
    local job_name=$2
    local memory=$3
    local cpus=$4
    local log_file="$LOG_DIR/${job_name}.log"
    local completion_file="$OUT_DIR/${job_name}.completed"

    echo "Submitting job for $script with name $job_name" | tee -a "$LOG_FILE"
    local job_id=$(bsub -cwd $(pwd) -R "rusage[mem=$memory]" -q $QUEUE_NAME -n "$cpus" -J "$job_name" -o "$log_file" "bash $script" 2>&1 | grep -oP '(?<=<)\d+(?=>)')
    if [ -z "$job_id" ]; then
        echo "Failed to submit job $job_name" | tee -a "$LOG_FILE"
        return 1
    fi

    echo "Job $job_name submitted with ID $job_id. Waiting for completion..." | tee -a "$LOG_FILE"
    # Wait for the job to finish
    bwait -w "done($job_id)"

    # Check the final status of the job
    local final_status=$(bjobs -noheader -o stat $job_id 2>/dev/null || echo "NOT FOUND")
    if [[ "$final_status" == "DONE" ]]; then
        echo "Job $job_name completed successfully." | tee -a "$LOG_FILE"
        # Create a completion file to indicate success
        touch "$completion_file"
        return 0
    else
        echo "Job $job_name failed or exited with status $final_status." | tee -a "$LOG_FILE"
        return 1
    fi
}



# Function to check if a step is completed
check_completed() {
    if [ -f "$OUT_DIR/$1.completed" ]; then
        echo "$1 already completed, skipping." | tee -a "$LOG_FILE"
        return 0
    else
        return 1
    fi
}

# Function to mark a step as completed
mark_completed() {
    touch "$OUT_DIR/$1.completed"
}

# Execute each step in the pipeline
execute_step() {
    local step_script=$1
    local step_name=$2
    local step_mem=$3
    local step_cpus=$4

    if ! check_completed "$step_name"; then
        echo "Executing $step_name..." | tee -a "$LOG_FILE"
        if submit_and_wait "$step_script" "${step_name}_job" "$step_mem" "$step_cpus"; then
            echo "$step_name completed successfully." | tee -a "$LOG_FILE"
            mark_completed "$step_name"
        else
            echo "Error executing $step_name. Check logs for details." | tee -a "$LOG_FILE"
            exit 1
        fi
    else
        echo "Skipping $step_name step as it is already completed."
    fi
}

# Step 1: Filter VCF files
execute_step "filter_vcf.sh" "filter_vcf" "$FILTER_JOB_MEM" "$FILTER_JOB_CPUS"

# Step 2: Merge VCF files
execute_step "merge_vcf.sh" "merge_vcf" "$MERGE_JOB_MEM" "$MERGE_JOB_CPUS"

# ... similar process for cellsnp.sh, vireo.sh, and souporcell.sh ...

# Example for cellsnp.sh
execute_step "cellsnp_job.sh" "cellsnp" "$CELLSNP_JOB_MEM" "$CELLSNP_JOB_CPUS"

# Example for vireo.sh
execute_step "vireo_job.sh" "vireo" "$VIREO_JOB_MEM" "$VIREO_JOB_CPUS"

# Example for souporcell.sh
# execute_step "souporcell.sh" "souporcell" "$SOUPORCELL_JOB_MEM" "$SOUPORCELL_JOB_CPUS"

echo "Pipeline execution completed." | tee -a "$LOG_FILE"
