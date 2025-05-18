#!/bin/bash

#PBS -j "job_array[1-4]%100"                   # Job array task range with a maximum of 5 concurrent tasks
#PBS -oo output.%I.txt                         # Output file pattern with task ID
#PBS -eo error.%I.txt                          # Error file pattern with task ID
#PBS -n 10                                      # Number of CPU cores required per task
#PBS -R "span[hosts=1]"                        # Number of hosts required per task
#PBS -R "rusage[mem=100G]"                   # Memory requirement per task (100G = 100,000MB)
#PBS -q verylong                               # Queue selection (verylong in this case)

# Set the number of task IDs to submit together
BATCH_SIZE=4

file='/omics/odcf/analysis/OE0536_projects/heroes-aya/scRNA-Seq/Demultiplexing_results/vireo/input_files/samples.txt' 

# Read the sample names and directories from sample.txt file
while IFS=$'\t' read -r dir  sample; do
  
  echo $sample
  # here we run vireoSNP tool to demultiplex 10x scRNAseq samples without genotyping data
  OUT_DIR=/omics/odcf/analysis/OE0536_projects/heroes-aya/scRNA-Seq/Demultiplexing_results/vireo/results/${sample}/
  CELL=/omics/odcf/analysis/OE0536_projects/heroes-aya/scRNA-Seq/Demultiplexing_results/cellsnp/results/${sample}/cellSNP.cells.vcf.gz
  CELL_DATA=/omics/odcf/analysis/OE0536_projects/heroes-aya/scRNA-Seq/Demultiplexing_results/cellsnp/results/${sample}/
  GT=/home/i439h/projects/heroes-aya/scRNA-Seq/Demultiplexing_results/liftover/vcf_files/Pool1.vcf.gz 
  n_donor=4

  ## create results subfolder for the sample 
  if [ ! -d "$OUT_DIR" ]; then
    mkdir "$OUT_DIR"
  fi

  cd $OUT_DIR

  # load modules
  module load bcftools/1.9
  module load vcftools/0.1.16


  # filter donor.vcf.gz files for those maintained in cellsnp output
  bcftools view $GT -R $CELL -Oz -o input.vcf.gz


  #vireo -c $CELL_DATA -N $n_donor -o $OUT_DIR -p 20 --callAmbientRNAs  --randSeed=12
  vireo -c $CELL_DATA -d input.vcf.gz -o $OUT_DIR -N $n_donor --randSee=12 

    # Check if the current task ID is divisible by the batch size
    if (( LSB_JOBINDEX % BATCH_SIZE == 0 )); then
      # Wait for all the previously submitted tasks to finish
      wait
    fi

done < $file
