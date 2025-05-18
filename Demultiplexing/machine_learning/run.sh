#!/bin/bash

vcf=/home/i439h/projects/pools/AG_Thongjuea/Result/10x_multiomics/Demultiplexing_results/vcf_files/pool3.merged.vcf
bam=/home/i439h/projects/pools/AG_Thongjuea/Result/10x_multiomics/CellRanger_ARC/pool3/outs/gex_possorted_bam.bam
barcode=/home/i439h/projects/pools/AG_Thongjuea/Result/10x_multiomics/CellRanger_ARC/pool3/outs/filtered_feature_bc_matrix/barcodes.tsv.gz
output=/home/i439h/projects/pools/AG_Thongjuea/Code/10x_multiomics/Sadeghi/Demultiplexing/machine_learning/test_output.csv

module load python/3.10.1 
source /home/i439h/projects/pools/AG_Thongjuea/Software/10x_multiomics/python/pyth_env/bin/activate
python demultiplex.py --vcf_file $vcf --bam_file $bam --barcode_file $barcode --n_clusters 4 --output_file $output
