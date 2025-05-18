#!/bin/bash
## demultiplex single cell multiomics-RNA data with unique genotyping data
# load modules
module load anaconda3/2021.05
module load singularity/3.10.5

# Get the pool and ID from command-line arguments
sample_number=$1
ID=$2

sample="pool${sample_number}_$ID"
pool="pool${sample_number}"

OUT_DIR=/results/multiomix/RNA/WGS/unique_snps/${sample}
BAM=/CellRanger_ARC/${ID}_${pool}_hg38/outs/gex_possorted_bam.bam
BARCODE=/CellRanger_ARC/${ID}_${pool}_hg38/outs/filtered_feature_bc_matrix/barcodes.tsv.gz
FASTA=/fasta/genome.fa
GT=/vcf_files/$pool.unique2.vcf

mkdir -p /omics/odcf/analysis/OE0290_projects/heroes-aya_pools/temp_analysis/scRNA-Seq/Demultiplexing_results/souporcell/results/multiomix/RNA/WGS/unique_snps/${sample}
cd /omics/odcf/analysis/OE0290_projects/heroes-aya_pools/temp_analysis/scRNA-Seq/Demultiplexing_results/souporcell/codes/singularity

singularity exec --bind /omics/odcf/analysis/OE0290_projects/heroes-aya_pools/temp_analysis/scRNA-Seq/Demultiplexing_results/souporcell/results:/results \
    --bind  /omics/odcf/analysis/OE0290_projects/heroes-aya_pools/AG_Thongjuea/Result/10x_multiomics/CellRanger_ARC:/CellRanger_ARC \
    --bind /home/i439h/projects/heroes-aya-pools/temp_analysis/scRNA-Seq/Demultiplexing_results/souporcell/input_files/vcf_files:/vcf_files \
    --bind /omics/odcf/analysis/OE0290_projects/heroes-aya_pools/temp_analysis/scRNA-Seq/CellRanger_v7.1_results/References/10x/refdata-gex-GRCh38-2020-A/fasta:/fasta souporcell_latest.sif souporcell_pipeline.py \
    -i $BAM \
    -f $FASTA \
    -t 20 \
    -o $OUT_DIR \
    -b $BARCODE \
    -k 4 \
    --known_genotypes $GT 
   ## --known_genotypes_sample_names I005_014 I010_013 I020_001 I023_037 
    
