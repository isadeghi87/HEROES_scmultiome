#!/bin/bash

# Define an array of IDs
declare -a IDs=(pool3)


vireo="/home/i439h/projects/pools/AG_Thongjuea/Code/10x_multiomics/Sadeghi/Demultiplexing/demuxafy/merged/no_genotype/vireo.sh"
demux="/home/i439h/projects/pools/AG_Thongjuea/Code/10x_multiomics/Sadeghi/Demultiplexing/demuxafy/merged/no_genotype/demuxalot.sh"
souporcell="/home/i439h/projects/pools/AG_Thongjuea/Code/10x_multiomics/Sadeghi/Demultiplexing/demuxafy/merged/no_genotype/souporcell.sh"

# Loop over each pool and ID, submitting a job for each
#for i in {1}; do
    #bsub -R "rusage[mem=100G]" -q long -n 10 -J "demux${i}" "$demux $i ${IDs[$i-1]}"
    bsub -R "rusage[mem=100G]" -q verylong -n 20 -J "vireo" "$vireo pool3"
    bsub -R "rusage[mem=100G]" -q verylong -n 20 -J "soup" "$souporcell pool3"
#done
