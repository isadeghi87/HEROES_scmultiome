#!/bin/bash

#BSUB -J cellranger_array[1-2]   # Job name and array range
#BSUB -o cellranger_output_%I.out   # Output file
#BSUB -e cellranger_error_%I.err   # Error file
#BSUB -n 30                        # Number of cores
#BSUB -R "rusage[mem=150G]" 
#BSUB -q verylong

# Load the Cell Ranger module
module add cellranger-arc/2.0.2 

# Set the path to the reference package
reference="/home/i439h/projects/pools/AG_Thongjuea/Dataset/10x_multiomics/references/human/refdata-cellranger-arc-GRCh38-2020-A-2.0.0"

# Define an array of samples
samples=(pool{23..23})

# Set the path to the output directory
sample=${samples[$((LSB_JOBINDEX-1))]}
lib_file="/home/i439h/projects/pools/AG_Thongjuea/Code/10x_multiomics/Sadeghi/cellRanger_arc/${sample}.csv"
results="/home/i439h/projects/pools/AG_Thongjuea/Result/10x_multiomics/CellRanger_ARC/"

# Use the LSF job array index to get the right sample
echo $sample
echo $lib_file

# Check if the BAM files already exist
bam_gex="${results}${sample}/outs/gex_possorted_bam.bam"
bam_atac="${results}${sample}/outs/atac_possorted_bam.bam"

if [ -f "$bam_gex" ] && [ -f "$bam_atac" ]; then
  echo "BAM files for $sample already exist. Skipping Cell Ranger."
else
  # Run cellranger arc for the current sample
  cellranger-arc count --id="$sample" \
    --reference="$reference" \
    --libraries=$lib_file \
    --localcores=20

  mv $sample $results
fi

# Submit the demultiplexing job for this pool
demux=/home/i439h/projects/pools/AG_Thongjuea/Code/10x_multiomics/Sadeghi/cellranger_demux_pipeline/combined_demux.sh
bsub -w "done($LSB_JOBID)" -R "rusage[mem=150G]" -q verylong -n 30 -J "demux_$sample" $demux $sample
