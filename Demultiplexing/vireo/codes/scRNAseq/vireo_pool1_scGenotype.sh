#!/bin/bash

module load bcftools/1.9
module load vcftools/0.1.16

# here we run vireoSNP tool to demultiplex 10x scRNAseq samples without genotyping data
sample='pool1_1016711'
OUT_DIR=/omics/odcf/analysis/OE0536_projects/heroes-aya/scRNA-Seq/Demultiplexing_results/vireo/results/${sample}/scGenotype
CELL=/omics/odcf/analysis/OE0536_projects/heroes-aya/scRNA-Seq/Demultiplexing_results/cellsnp/results/${sample}/scGenotype/cellSNP.cells.vcf.gz
CELL_DATA=/omics/odcf/analysis/OE0536_projects/heroes-aya/scRNA-Seq/Demultiplexing_results/cellsnp/results/${sample}/scGenotype/
GT=/omics/odcf/analysis/OE0536_projects/heroes-aya/scWGS_PTA/calledVariants/freebayes/scWGS_pool1.vcf
n_donor=4

## create results subfolder for the sample 
if [ ! -d "$OUT_DIR" ]; then
  mkdir "$OUT_DIR"
fi

cd $OUT_DIR

# filter donor.vcf.gz files for those maintained in cellsnp output


#vireo -c $CELL_DATA -N $n_donor -o $OUT_DIR -p 20 --callAmbientRNAs  --randSeed=12
vireo -c $CELL_DATA -d $GT -o $OUT_DIR -N $n_donor --randSee=12 --genoTag=GT 