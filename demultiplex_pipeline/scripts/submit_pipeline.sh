#!/bin/bash

# Load any necessary modules or environment variables here
# module load your_module
source  /home/i439h/projects/heroes-aya-pools/AG_Thongjuea/Code/10x_multiomics/Sadeghi/demultiplex_pipeline/scripts/config.sh

# Define job submission parameters
PIPELINE_SCRIPT_PATH="/home/i439h/projects/heroes-aya-pools/AG_Thongjuea/Code/10x_multiomics/Sadeghi/demultiplex_pipeline/scripts/run_pipeline.sh"
OUTPUT_FILE="$LOG_DIR/${JOB_NAME}_output.log"
LOG_DIR="/home/i439h/projects/heroes-aya-pools/AG_Thongjuea/Code/10x_multiomics/Sadeghi/demultiplex_pipeline/output/logs"
ERROR_FILE="$LOG_DIR/${JOB_NAME}_error.log"

QUEUE_NAME="long"  # replace with your LSF queue name
JOB_NAME="demultiplex_pipeline"
JOB_MEMORY="100G"  # Adjust as needed for each job
JOB_CPUS=10
working_dir=$(pwd)
# Make sure the log directory exists
if [ ! -d "$LOG_DIR" ]; then
    mkdir -p "$LOG_DIR"
fi

# Submit the pipeline script as a job
bsub -cwd "$working_dir" -R "rusage[mem=$JOB_MEMORY]" -q "$QUEUE_NAME" -n "$JOB_CPUS" -J "pipeline_job" -o "$LOG_FILE" "bash $PIPELINE_SCRIPT_PATH"