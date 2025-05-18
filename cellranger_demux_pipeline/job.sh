#!/bin/bash

# Define an array of samples
samples=(pool{21..22})
cellranger=/home/i439h/projects/pools/AG_Thongjuea/Code/10x_multiomics/Sadeghi/cellranger_demux_pipeline/cellranger_script.sh
# Loop over each sample and submit the Cell Ranger job
for i in "${!samples[@]}"; do
  sample="${samples[$i]}"
  bsub -J "cellranger_${i}" -o "cellranger_output_${i}.out" -e "cellranger_error_${i}.err" -n 30 -R "rusage[mem=150G]" -q verylong $cellranger
done
