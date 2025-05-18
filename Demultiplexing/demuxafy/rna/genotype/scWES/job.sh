#!/bin/bash

# Define an array of IDs
declare -a IDs=('pool3' 'pool4')

vireo="/home/i439h/projects/pools/AG_Thongjuea/Code/10x_multiomics/Sadeghi/Demultiplexing/demuxafy/rna/genotype/scWES/vireo.sh"
#demux="/home/i439h/projects/pools/AG_Thongjuea/Code/10x_multiomics/Sadeghi/Demultiplexing/demuxafy/rna/genotype/scWES/demuxalot.sh"
souporcell="/home/i439h/projects/pools/AG_Thongjuea/Code/10x_multiomics/Sadeghi/Demultiplexing/demuxafy/rna/genotype/scWES/souporcell.sh"

# Loop over each pool and ID, submitting a job for each
for i in {1..2}; do
    #bsub -R "rusage[mem=100G]" -q long -n 10 -J "demux${i}" "$demux $i ${IDs[$i-1]}"
    bsub -R "rusage[mem=100G]" -q long -n 10 -J "vireo${i}" "$vireo ${IDs[$i-1]}"
    #bsub -R "rusage[mem=100G]" -q long -n 10 -J "soup${i}" "$souporcell ${IDs[$i-1]}"
done
