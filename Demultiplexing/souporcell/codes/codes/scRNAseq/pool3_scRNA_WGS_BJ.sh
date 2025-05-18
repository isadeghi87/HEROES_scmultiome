#!/bin/bash
## demultiplex single cell-RNA data using bulkWGS snps from basejumper (BJ)
# load modules
module load anaconda3/2021.05
module load singularity/3.10.5

ID='1016715'
pool='pool3'
sample=${pool}_${ID}
OUT_DIR=/results/scRNAseq/${sample}/WGS_BJ
BAM=/CellRanger_Results/${ID}_${pool}_hg38/outs/possorted_genome_bam.bam
BARCODE=/CellRanger_Results/${ID}_${pool}_hg38/outs/filtered_feature_bc_matrix/barcodes.tsv.gz
FASTA=/fasta/genome.fa
GT=/vcf_files/pool4.qual.vcf



singularity exec --bind /omics/odcf/analysis/OE0290_projects/heroes-aya_pools/temp_analysis/scRNA-Seq/Demultiplexing_results/souporcell/results:/results \
    --bind /omics/odcf/analysis/OE0290_projects/heroes-aya_pools/temp_analysis/scRNA-Seq/CellRanger_v7.1_results/CellRanger_Results:/CellRanger_Results \
    --bind /omics/odcf/analysis/OE0290_projects/heroes-aya_pools/temp_analysis/scRNA-Seq/Demultiplexing_results/souporcell/input_files/vcf_files:/vcf_files \
    --bind /omics/odcf/analysis/OE0290_projects/heroes-aya_pools/temp_analysis/scRNA-Seq/CellRanger_v7.1_results/References/10x/refdata-gex-GRCh38-2020-A/fasta:/fasta souporcell_latest.sif souporcell_pipeline.py \
      -i $BAM \
    -f $FASTA \
    -t 20 \
    -o $OUT_DIR \
    -b $BARCODE \
    -k 4 \
    --known_genotypes $GT 
