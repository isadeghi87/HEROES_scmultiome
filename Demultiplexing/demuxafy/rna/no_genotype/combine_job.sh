#!/bin/bash

# Define an array of IDs
declare -a IDs=(pool{1..22})


combine=/home/i439h/projects/pools/AG_Thongjuea/Code/10x_multiomics/Sadeghi/Demultiplexing/demuxafy/rna/no_genotype/combine_results.sh

# Loop over each pool and ID, submitting a job for each
for i in {1..22}; do
    bsub -R "rusage[mem=100G]" -q long -n 10 -J "comb${i}" "$combine  ${IDs[$i-1]}"
done
