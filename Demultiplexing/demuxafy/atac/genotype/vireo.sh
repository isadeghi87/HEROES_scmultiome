#!/bin/bash

module load anaconda3/2021.05
module load singularity/3.10.5
module load bcftools/1.9

# Get the pool and ID from command-line arguments
pool="${1:?Pool not specified}"
tool='vireo'
OUT_DIR=/$tool/$pool
BAM=/CellRanger_ARC/${pool}/outs/atac_possorted_bam.bam
BARCODES=/CellRanger_ARC/${pool}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz
#VCF=/vcf_files/${pool}.merged.vcf.gz
VCF=/vcf_files/${pool}.merged.vcf.gz
FIELD='GT' 

result=/omics/odcf/analysis/OE0290_projects/heroes-aya_pools/AG_Thongjuea/Result/10x_multiomics/Demultiplexing_results/demuxafy/ATAC/genotype/$tool/$pool

cd /home/i439h/projects/pools/AG_Thongjuea/Software/10x_multiomics/demuxafy

# create results subfolder for the sample 
if [ ! -d "$result" ]; then  
  mkdir -p "$result"
fi

## run cellsnp-lite
singularity exec --bind /omics/odcf/analysis/OE0290_projects/heroes-aya_pools/AG_Thongjuea/Result/10x_multiomics/Demultiplexing_results/demuxafy/ATAC/genotype/$tool:/$tool \
    --bind /home/i439h/projects/pools/AG_Thongjuea/Result/10x_multiomics/CellRanger_ARC:/CellRanger_ARC \
    --bind /home/i439h/projects/pools/AG_Thongjuea/Code/10x_multiomics/Sadeghi/Demultiplexing/liftover/vcf_files:/vcf_files \
    --bind /home/i439h/projects/pools/AG_Thongjuea/Code/10x_multiomics/Sadeghi/Demultiplexing/demuxafy/atac:/atac Demuxafy.sif cellsnp-lite  \
        -b $BARCODES \
        -R $VCF \
        -O $OUT_DIR \
        --cellTAG CB \
        --UMItag None \
        -s $BAM \
        -p 20 \
        --minMAF 0.1 \
        --minCOUNT 20 \
        --gzip \
        --genotype

VCF=/home/i439h/projects/pools/AG_Thongjuea/Code/10x_multiomics/Sadeghi/Demultiplexing/liftover/vcf_files/${pool}.merged.vcf.gz
CELL=/omics/odcf/analysis/OE0290_projects/heroes-aya_pools/AG_Thongjuea/Result/10x_multiomics/Demultiplexing_results/demuxafy/ATAC/genotype/$tool/$pool/cellSNP.cells.vcf.gz
out=/omics/odcf/analysis/OE0290_projects/heroes-aya_pools/AG_Thongjuea/Result/10x_multiomics/Demultiplexing_results/demuxafy/ATAC/genotype/$tool/$pool/
bcftools view $VCF -R $CELL -Oz -o ${out}/input.vcf.gz

singularity exec --bind /omics/odcf/analysis/OE0290_projects/heroes-aya_pools/AG_Thongjuea/Result/10x_multiomics/Demultiplexing_results/demuxafy/ATAC/genotype/$tool:/$tool \
    --bind /home/i439h/projects/pools/AG_Thongjuea/Result/10x_multiomics/Demultiplexing_results/vcf_files:/vcf_files \
    --bind /home/i439h/projects/pools/AG_Thongjuea/Code/10x_multiomics/Sadeghi/Demultiplexing/demuxafy/atac:/atac Demuxafy.sif vireo \
    -c $OUT_DIR \
    -d ${OUT_DIR}/input.vcf.gz \
    -o $OUT_DIR \
    --genoTag=GT \
    -p 20 \
    -N 4 \
    --forceLearnGT