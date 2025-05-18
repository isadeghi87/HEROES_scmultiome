#!/bin/bash

module load anaconda3/2021.05
source activate scavengers
cd /home/i439h/projects/heroes-aya/scRNA-Seq/Demultiplexing_results/souporcell/results/pool3_1016715
/home/i439h/projects/heroes-aya/scRNA-Seq/Demultiplexing_results/scAVENGERS/codes/scAVENGERS-0.2.0/scAVENGERS cluster -a alt.mtx -r ref.mtx -b barcodes.tsv -o ./ -k 4
/home/i439h/projects/heroes-aya/scRNA-Seq/Demultiplexing_results/scAVENGERS/codes/scAVENGERS-0.2.0/scAVENGERS assign -g gt_matrix.npz -i variant_index.npz -v variants.vcf.gz -r donor_variants.vcf.gz