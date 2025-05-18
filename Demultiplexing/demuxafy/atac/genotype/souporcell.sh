#!/bin/bash

module load anaconda3/2021.05
module load singularity/3.10.5
module load bcftools/1.9
module load R/4.3.0

# Get the pool and ID from command-line arguments
pool_number=$1
ID=$2
sample="pool${pool_number}_$ID"
pool="pool${pool_number}"
tool='souporcell'
OUT_DIR=/$tool/$pool
BAM=/CellRanger_ARC/${ID}_${pool}_hg38/outs/atac_possorted_bam.bam
BARCODES=/CellRanger_ARC/${ID}_${pool}_hg38/outs/filtered_feature_bc_matrix/barcodes.tsv.gz
VCF=/vcf_files/${pool}.merged.vcf
FASTA=/hg38/genome.fa
FIELD='GT' 
INDS=/atac/${pool}_donor_list.txt
N=4
THREADS=20
DIR="/home/i439h/projects/heroes-aya-pools/AG_Thongjuea"

# create results subfolder for each pool 
result=$DIR/Result/10x_multiomics/Demultiplexing_results/demuxafy/ATAC/genotype/$tool/$pool
if [ ! -d "$result" ]; then  
  mkdir -p "$result"
fi

singularity exec --bind $DIR/Result/10x_multiomics/Demultiplexing_results/demuxafy/ATAC/genotype/$tool:/$tool \
    --bind $DIR/Result/10x_multiomics/CellRanger_ARC/:/CellRanger_ARC \
    --bind $DIR/Result/10x_multiomics/Demultiplexing_results/vcf_files:/vcf_files \
    --bind $DIR/Code/10x_multiomics/Sadeghi/Demultiplexing/demuxafy/atac:/atac \
    --bind $DIR/Code/10x_multiomics/Sadeghi/Demultiplexing/references/hg38:/hg38 Demuxafy.sif souporcell_pipeline.py \
    -i $BAM \
    -b $BARCODES \
    -f $FASTA \
    -t $THREADS \
    -o $OUT_DIR \
    -k $N \
    --known_genotypes $VCF \
    --no_umi True \
    --skip_remap True
       
#Souporcell Summary
singularity exec --bind $DIR/Result/10x_multiomics/Demultiplexing_results/demuxafy/ATAC/genotype/$tool:/$tool Demuxafy.sif bash ${tool}_summary.sh $OUT_DIR/clusters.tsv > $DIR/Result/10x_multiomics/Demultiplexing_results/demuxafy/ATAC/genotype/$tool/$pool/${tool}_summary.tsv

#Correlating Cluster to Donor Reference SNP Genotypes
singularity exec --bind $DIR/Result/10x_multiomics/Demultiplexing_results/vcf_files:/vcf_files \
        --bind $DIR/Result/10x_multiomics/Demultiplexing_results/demuxafy/ATAC/genotype/$tool:/$tool Demuxafy.sif Assign_Indiv_by_Geno.R \
        -r $VCF -c $OUT_DIR/cluster_genotypes.vcf -o $OUT_DIR

echo "Pool: $pool"

Rscript souporcell_assign.R "$pool"