#!/bin/bash

# Define an array of IDs
declare -a IDs=("pool14")

vireo="/home/i439h/projects/pools/AG_Thongjuea/Code/10x_multiomics/Sadeghi/Demultiplexing/demuxafy/atac/genotype/vireo.sh"
demux="/home/i439h/projects/pools/AG_Thongjuea/Code/10x_multiomics/Sadeghi/Demultiplexing/demuxafy/atac/genotype/demuxalot.sh"
souporcell="/home/i439h/projects/pools/AG_Thongjuea/Code/10x_multiomics/Sadeghi/Demultiplexing/demuxafy/atac/genotype/souporcell.sh"

# Loop over each pool and ID, submitting a job for each
for i in {1..2}; do
    #bsub -R "rusage[mem=100G]" -q long -n 10 -J "demux${i}" "$demux $i ${IDs[$i-1]}"
    bsub -R "rusage[mem=150G]" -q verylong -n 20 -J "vireo${i}" "$vireo ${IDs[$i-1]}"
    #bsub -R "rusage[mem=100G]" -q long -n 30 -J "soup${i}" "$souporcell $i ${IDs[$i-1]}"
done
