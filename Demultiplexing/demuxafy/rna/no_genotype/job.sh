#!/bin/bash

# Define an array of IDs
declare -a IDs=(pool{9..22})


vireo="/home/i439h/projects/pools/AG_Thongjuea/Code/10x_multiomics/Sadeghi/Demultiplexing/demuxafy/rna/no_genotype/vireo.sh"
souporcell="/home/i439h/projects/pools/AG_Thongjuea/Code/10x_multiomics/Sadeghi/Demultiplexing/demuxafy/rna/no_genotype/souporcell.sh"

# Loop over each pool and ID, submitting a job for each
for i in {1..14}; do
    #bsub -R "rusage[mem=150G]" -q verylong -n 20 -J "vireo${i}" "$vireo  ${IDs[$i-1]}"
    bsub -R "rusage[mem=150G]" -q verylong -n 20 -J "soup${i}" "$souporcell  ${IDs[$i-1]}"
done
