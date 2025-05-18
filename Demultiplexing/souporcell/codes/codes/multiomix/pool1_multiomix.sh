#!/bin/bash
# load modules
module load anaconda3/2021.05


module load singularity/3.10.5
# here we run souporcell tool to demultiplex 10x multiomix samples without genotyping data
ID='1016711'
pool='pool1'
sample=${pool}_${ID}
OUT_DIR=/results/${sample}/multiomix/scRNAseq
BAM=/CellRanger_arc/${ID}_${pool}_hg38/outs/gex_possorted_bam.bam
BARCODE=/CellRanger_arc/${ID}_${pool}_hg38/outs/filtered_feature_bc_matrix/barcodes.tsv.gz
FASTA=/fasta/genome.fa
#GT=/freebayes/scWGS_pool1.vcf

## create results subfolder for the sample 
if [ ! -d "$OUT_DIR" ]; then
  mkdir -p "$OUT_DIR"
fi


singularity exec --bind /omics/odcf/analysis/OE0536_projects/heroes-aya/scRNA-Seq/Demultiplexing_results/souporcell/results:/results,/omics/odcf/analysis/OE0536_projects/heroes-aya/scMulti-omics/CellRanger_arc:/CellRanger_arc,/omics/odcf/analysis/OE0536_projects/heroes-aya/scRNA-Seq/CellRanger_v7.1_results/References/10x/refdata-gex-GRCh38-2020-A/fasta:/fasta souporcell_latest.sif souporcell_pipeline.py \
    -i $BAM \
    -f $FASTA \
    -t 20 \
    -o $OUT_DIR \
    -b $BARCODE \
    -k 4

