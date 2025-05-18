#!/bin/bash
## running demuxlet from demuxafy for pool3

module load anaconda3/2021.05
module load singularity/3.10.5

OUT_DIR=/demuxlet/pool3
BAM=/CellRanger_ARC/1016715_pool3_hg38/outs/gex_possorted_bam.bam
BARCODE=/CellRanger_ARC/1016715_pool3_hg38/outs/filtered_feature_bc_matrix/barcodes.tsv.gz
FIELD='GT' ## this might also be GT depending on the fields in your vcf
VCF=/vcf_files/pool3.merged.vcf
INDS=/demuxafy/donor_list.txt



singularity exec --bind /omics/odcf/analysis/OE0290_projects/heroes-aya_pools/AG_Thongjuea/Result/10x_multiomics/Demultiplexing_results/demuxafy/demuxlet:/demuxlet \
    --bind /home/i439h/projects/heroes-aya-pools/AG_Thongjuea/Result/10x_multiomics/CellRanger_ARC/:/CellRanger_ARC \
    --bind /home/i439h/projects/heroes-aya-pools/AG_Thongjuea/Result/10x_multiomics/Demultiplexing_results/vcf_files:/vcf_files \
    --bind /home/i439h/projects/heroes-aya-pools/AG_Thongjuea/Code/10x_multiomics/Sadeghi/Demultiplexing/demuxafy:/demuxafy Demuxafy.sif popscle dsc-pileup  \
        --sam $BAM \
        --vcf $VCF \
        --group-list $BARCODE \
        --out $OUT_DIR \
        --sm-list $INDS

# singularity exec --bind /omics/odcf/analysis/OE0290_projects/heroes-aya_pools/AG_Thongjuea/Result/10x_multiomics/Demultiplexing_results/demuxafy/demuxlet:/demuxlet \
#     --bind /home/i439h/projects/heroes-aya-pools/AG_Thongjuea/Result/10x_multiomics/CellRanger_ARC/:/CellRanger_ARC \
#     --bind /home/i439h/projects/heroes-aya-pools/AG_Thongjuea/Result/10x_multiomics/Demultiplexing_results/vcf_files:/vcf_files \
#     --bind /home/i439h/projects/heroes-aya-pools/AG_Thongjuea/Code/10x_multiomics/Sadeghi/Demultiplexing/demuxafy:/demuxafy \
#     --bind /home/i439h/projects/heroes-aya-pools/AG_Thongjuea/Code/10x_multiomics/Sadeghi/Demultiplexing/references/hg38:/hg38 Demuxafy.sif popscle demuxlet --plp $OUT_DIR \
#         --vcf $VCF --field $FIELD --group-list $BARCODES --geno-error-coeff 1.0 --geno-error-offset 0.05 --out $OUT_DIR --sm-list $INDS


# singularity exec --bind /omics/odcf/analysis/OE0290_projects/heroes-aya_pools/AG_Thongjuea/Result/10x_multiomics/Demultiplexing_results/demuxafy/demuxlet:/demuxlet Demuxafy.sif bash demuxlet_summary.sh $OUT_DIR/demuxlet.best > /omics/odcf/analysis/OE0290_projects/heroes-aya_pools/AG_Thongjuea/Result/10x_multiomics/Demultiplexing_results/demuxafy/demuxlet/pool3/demuxlet_summary.tsv
