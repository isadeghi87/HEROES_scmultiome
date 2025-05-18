#!/bin/bash

# Define an array of IDs
declare -a IDs=(pool{9..22})

# Path to the combined script
combined_script="/home/i439h/projects/pools/AG_Thongjuea/Code/10x_multiomics/Sadeghi/Demultiplexing/demuxafy/automated/combined_demuxafy.sh"

# Loop over each pool and submit a job for each
for i in {0..13}; do
    bsub -R "rusage[mem=150G]" -q verylong -n 30 -J "demux_${IDs[$i]}" "$combined_script ${IDs[$i]}"
done
