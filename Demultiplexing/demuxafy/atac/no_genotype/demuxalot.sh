#!/bin/bash

module load anaconda3/2021.05
module load singularity/3.10.5
module load bcftools/1.9

# Get the pool and ID from command-line arguments
pool_number=$1
ID=$2
sample="pool${pool_number}_$ID"
pool="pool${pool_number}"
tool='demuxalot'
OUT_DIR=/$tool/$pool

BAM=/CellRanger_ARC/${ID}_${pool}_hg38/outs/gex_possorted_bam.bam
BARCODES=/CellRanger_ARC/${ID}_${pool}_hg38/outs/filtered_feature_bc_matrix/barcodes.tsv.gz
VCF=/vcf_files/${pool}.merged.vcf
FASTA=/hg38/genome.fa

#Get sample names and save to a file
# vcf=/home/i439h/projects/heroes-aya-pools/AG_Thongjuea/Result/10x_multiomics/Demultiplexing_results/vcf_files/${pool}.merged.vcf
# bcftools query -l "$vcf" > /omics/odcf/analysis/OE0290_projects/heroes-aya_pools/AG_Thongjuea/Result/10x_multiomics/Demultiplexing_results/demuxafy/ATAC/no_genotype/${pool}_donor_list.txt

INDS=/atac/${pool}_donor_list.txt

# create results subfolder for each pool 
result=/omics/odcf/analysis/OE0290_projects/heroes-aya_pools/AG_Thongjuea/Result/10x_multiomics/Demultiplexing_results/demuxafy/ATAC/no_genotype/$tool/$pool
if [ ! -d "$result" ]; then  
  mkdir -p "$result"
fi


singularity exec --bind /omics/odcf/analysis/OE0290_projects/heroes-aya_pools/AG_Thongjuea/Result/10x_multiomics/Demultiplexing_results/demuxafy/ATAC/no_genotype/$tool:/$tool \
    --bind /home/i439h/projects/heroes-aya-pools/AG_Thongjuea/Result/10x_multiomics/CellRanger_ARC/:/CellRanger_ARC \
    --bind /home/i439h/projects/heroes-aya-pools/AG_Thongjuea/Result/10x_multiomics/Demultiplexing_results/vcf_files:/vcf_files \
    --bind /home/i439h/projects/heroes-aya-pools/AG_Thongjuea/Code/10x_multiomics/Sadeghi/Demultiplexing/demuxafy/atac:/atac \
    --bind /home/i439h/projects/heroes-aya-pools/AG_Thongjuea/Code/10x_multiomics/Sadeghi/Demultiplexing/references/hg38:/hg38 Demuxafy.sif Demuxalot.py  \
        -b $BARCODES \
        -a $BAM \
        -n $INDS \
        -v $VCF \
        -o $OUT_DIR \
        -r True


singularity exec --bind /omics/odcf/analysis/OE0290_projects/heroes-aya_pools/AG_Thongjuea/Result/10x_multiomics/Demultiplexing_results/demuxafy/ATAC/no_genotype/$tool:/$tool Demuxafy.sif bash ${tool}_summary.sh $OUT_DIR/assignments_refined.tsv.gz > /omics/odcf/analysis/OE0290_projects/heroes-aya_pools/AG_Thongjuea/Result/10x_multiomics/Demultiplexing_results/demuxafy/ATAC/no_genotype/$tool/$pool/${tool}_summary.tsv
