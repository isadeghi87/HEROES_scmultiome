#!/bin/bash
module load vartrix/1.1.22

vcf=/home/i439h/projects/heroes-aya-pools/AG_Thongjuea/Result/bulkWGS/mitochondria/vcfs/pool4.merged.vcf.gz
fasta="/home/i439h/projects/heroes-aya-pools/AG_Thongjuea/Code/10x_multiomics/Sadeghi/Demultiplexing/references/hg38/genome.fa"
bam=/home/i439h/projects/heroes-aya-pools/AG_Thongjuea/Result/10x_multiomics/CellRanger_ARC/1016717_pool4_hg38/outs/gex_possorted_bam.bam
outdir=/home/i439h/projects/heroes-aya-pools/AG_Thongjuea/Result/10x_multiomics/Demultiplexing_results/vartrix/mito_scRNA/pool4
barcode=/home/i439h/projects/heroes-aya-pools/AG_Thongjuea/Result/10x_multiomics/CellRanger_ARC/1016717_pool4_hg38/outs/filtered_feature_bc_matrix/barcodes.tsv.gz
mkdir -p $outdir

vartrix \
    --bam $bam \
    --fasta $fasta \
    --vcf $vcf \
    --out-matrix $outdir/matrix.mtx \
    --cell-barcodes $barcode \
    --threads 20 \
    --mapq 10 \
    --umi
