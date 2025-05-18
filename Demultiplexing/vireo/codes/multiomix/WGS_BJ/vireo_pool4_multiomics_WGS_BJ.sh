#!/bin/bash
source ~/.bash_profile
module load bcftools/1.9
module load vcftools/0.1.16

# here we run vireoSNP tool to demultiplex 10x scRNAseq samples without genotyping data
sample='pool4_1016717'
OUT_DIR=/home/i439h/projects/heroes-aya-pools/temp_analysis/scRNA-Seq/Demultiplexing_results/vireo/results/${sample}/multiomics/RNA/WGS_BJ
CELL=/home/i439h/projects/heroes-aya-pools/temp_analysis/scRNA-Seq/Demultiplexing_results/cellsnp/results/${sample}/WGS_BJ/cellSNP.cells.vcf.gz
CELL_DATA=/home/i439h/projects/heroes-aya-pools/temp_analysis/scRNA-Seq/Demultiplexing_results/cellsnp/results/${sample}/WGS_BJ
GT=/home/i439h/projects/heroes-aya-pools/temp_analysis/scRNA-Seq/Demultiplexing_results/souporcell/input_files/vcf_files/pool4.filtered,.vcf.gz

n_donor=4

## create results subfolder for the sample 
if [ ! -d "$OUT_DIR" ]; then
  mkdir -p "$OUT_DIR"
fi

cd $OUT_DIR

# filter donor.vcf.gz files for those maintained in cellsnp output
bcftools view $GT -R $CELL -Oz -o input.vcf.gz

#vireo -c $CELL_DATA -N $n_donor -o $OUT_DIR -p 20 --callAmbientRNAs  --randSeed=12
vireo -c $CELL_DATA -d input.vcf.gz -o $OUT_DIR -N $n_donor --randSee=12 --genoTag=GT