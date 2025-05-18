#!/bin/bash

# Define an array of IDs
declare -a IDs=("1016711" "1016713" "1016715" "1016717")

SCRIPT_PATH="/home/i439h/projects/heroes-aya-pools/temp_analysis/scRNA-Seq/Demultiplexing_results/cellsnp/codes/multiomix/1000genome/cellsnp_multomix_RNA_1000genome.sh"

# Loop over each pool and ID, submitting a job for each
for i in {1..4}; do
    bsub -R "rusage[mem=200G]" -q long -n 10 -J "cell${i}_1000genome" \
    "$SCRIPT_PATH $i ${IDs[$i-1]}"
done
