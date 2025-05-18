#!/bin/bash
## demultiplex single cell multiomics-RNA data using bulkWGS
# load modules
module load anaconda3/2021.05
module load singularity/3.10.5

ID='1016715'
pool='pool3'
sample=${pool}_${ID}
OUT_DIR=/results/multiomix/RNA/with_bulkWGS/${sample}
BAM=/CellRanger_arc/${ID}_${pool}_hg38/outs/gex_possorted_bam.bam
BARCODE=/CellRanger_arc/${ID}_${pool}_hg38/outs/filtered_feature_bc_matrix/barcodes.tsv.gz
FASTA=/fasta/genome.fa
GT=/WGS/pool3.vcf



singularity exec --bind /omics/odcf/analysis/OE0290_projects/heroes-aya_pools/temp_analysis/scRNA-Seq/Demultiplexing_results/souporcell/results:/results \
    --bind /omics/odcf/analysis/OE0290_projects/heroes-aya_pools/temp_analysis/scMulti-omics/CellRanger_arc:/CellRanger_arc \
    --bind /omics/odcf/analysis/OE0290_projects/heroes-aya_pools/temp_analysis/scRNA-Seq/Demultiplexing_results/liftover/vcf_files/WGS:/WGS \
    --bind /omics/odcf/analysis/OE0290_projects/heroes-aya_pools/temp_analysis/scRNA-Seq/CellRanger_v7.1_results/References/10x/refdata-gex-GRCh38-2020-A/fasta:/fasta souporcell_latest.sif souporcell_pipeline.py \
    -i $BAM \
    -f $FASTA \
    -t 20 \
    -o $OUT_DIR \
    -b $BARCODE \
    -k 4 \
    --known_genotypes $GT \
    --known_genotypes_sample_names I010_013 I005_014 I020_001 I023_037 \
    
