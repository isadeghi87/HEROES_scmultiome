#!/bin/bash

# Define an array of IDs
declare -a IDs=("1016711" "1016713" "1016715" "1016717")

SCRIPT_PATH="/home/i439h/projects/heroes-aya-pools/AG_Thongjuea/Code/10x_multiomics/Sadeghi/Demultiplexing/cellsnp/codes/multiomix/WGS_BJ/q2000_dp100/cellsnp_multomix_RNA_WGS_q2000_dp100.sh"

# Loop over each pool and ID, submitting a job for each
for i in {1..4}; do
    bsub -R "rusage[mem=200G]" -q long -n 10 -J "cell${i}_q2000" \
    "$SCRIPT_PATH $i ${IDs[$i-1]}"
done
