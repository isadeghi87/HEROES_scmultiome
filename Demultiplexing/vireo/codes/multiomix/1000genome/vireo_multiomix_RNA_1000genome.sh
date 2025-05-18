#!/bin/bash
source ~/.bash_profile
module load bcftools/1.9
module load vcftools/0.1.16

# here we run vireoSNP tool to demultiplex 10x scRNAseq samples with 1000genome genotyping data
# Get the pool and ID from command-line arguments
sample_number=$1
ID=$2

sample="pool${sample_number}_$ID"
pool="pool${sample_number}"
OUT_DIR=/home/i439h/projects/heroes-aya-pools/temp_analysis/scRNA-Seq/Demultiplexing_results/vireo/results/multiomix/RNA/1000genome/${sample}
CELL=/home/i439h/projects/heroes-aya-pools/temp_analysis/scRNA-Seq/Demultiplexing_results/cellsnp/results/multiomix/RNA/1000genome/${sample}/cellSNP.cells.vcf.gz
CELL_DATA=/home/i439h/projects/heroes-aya-pools/temp_analysis/scRNA-Seq/Demultiplexing_results/cellsnp/results/multiomix/RNA/1000genome/${sample}/
#GT=/home/i439h/projects/heroes-aya-pools/temp_analysis/scRNA-Seq/Demultiplexing_results/cellsnp/input_files/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz

n_donor=4

## create results subfolder for the sample 
if [ ! -d "$OUT_DIR" ]; then
  mkdir -p "$OUT_DIR"
fi

cd $OUT_DIR

# filter donor.vcf.gz files for those maintained in cellsnp output
#bcftools view $GT -R $CELL -Oz -o input.vcf.gz

#vireo -c $CELL_DATA -N $n_donor -o $OUT_DIR -p 20 --callAmbientRNAs  --randSeed=12
vireo -c $CELL_DATA -o $OUT_DIR -N $n_donor --randSee=12 -p 20