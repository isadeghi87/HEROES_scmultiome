#!/bin/bash
# load modules
module load anaconda3/2021.05


module load singularity/3.10.5
# here we run souporcell tool to demultiplex 10x scRNAseq samples without genotyping data
sample='pool1_1016711'
pool='Pool1'
OUT_DIR=/results/${sample}/scGenotype
BAM=/input_files/${sample}/genome_bam.bam
BARCODE=/input_files/${sample}/barcodes.tsv
FASTA=/fasta/genome.fa
GT=/freebayes/scWGS_pool1.vcf

## create results subfolder for the sample 
if [ ! -d "$OUT_DIR" ]; then
  mkdir -p "$OUT_DIR"
fi


singularity exec --bind /omics/odcf/analysis/OE0536_projects/heroes-aya/scRNA-Seq/Demultiplexing_results/souporcell/results:/results,/omics/odcf/analysis/OE0536_projects/heroes-aya/scRNA-Seq/Demultiplexing_results/souporcell/input_files:/input_files,/omics/odcf/analysis/OE0536_projects/heroes-aya/scRNA-Seq/CellRanger_v7.1_results/References/10x/refdata-gex-GRCh38-2020-A/fasta:/fasta,/omics/odcf/analysis/OE0536_projects/heroes-aya/scWGS_PTA/calledVariants/freebayes:/freebayes souporcell_latest.sif souporcell_pipeline.py \
    -i $BAM \
    -f $FASTA \
    -t 20 \
    -o $OUT_DIR \
    -b $BARCODE \
    --umi_tag UB \
    --cell_tag CB \
    --known_genotypes $GT \
    -k 4