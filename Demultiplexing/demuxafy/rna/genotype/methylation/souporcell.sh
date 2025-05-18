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
BAM=/CellRanger_ARC/${ID}_${pool}_hg38/outs/gex_possorted_bam.bam
BARCODES=/CellRanger_ARC/${ID}_${pool}_hg38/outs/filtered_feature_bc_matrix/barcodes.tsv.gz
#VCF=/vcf_files/${pool}.merged.vcf
VCF=/vcf_files/I070_032_tumor001-01_control001-01_hg38.vcf.gz
FASTA=/hg38/genome.fa
FIELD='GT' 
INDS=/rna/${pool}_donor_list.txt
N=4
THREADS=20


# create results subfolder for each pool 
result=/omics/odcf/analysis/OE0290_projects/heroes-aya_pools/AG_Thongjuea/Result/10x_multiomics/Demultiplexing_results/demuxafy/RNA/$tool/$pool
if [ ! -d "$result" ]; then  
  mkdir -p "$result"
fi

singularity exec --bind /omics/odcf/analysis/OE0290_projects/heroes-aya_pools/AG_Thongjuea/Result/10x_multiomics/Demultiplexing_results/demuxafy/RNA/$tool:/$tool \
    --bind /home/i439h/projects/heroes-aya-pools/AG_Thongjuea/Result/10x_multiomics/CellRanger_ARC/:/CellRanger_ARC \
    --bind /home/i439h/projects/heroes-aya-pools/AG_Thongjuea/Result/10x_multiomics/Demultiplexing_results/vcf_files:/vcf_files \
    --bind /home/i439h/projects/heroes-aya-pools/AG_Thongjuea/Code/10x_multiomics/Sadeghi/Demultiplexing/demuxafy/rna:/rna \
    --bind /home/i439h/projects/heroes-aya-pools/AG_Thongjuea/Code/10x_multiomics/Sadeghi/Demultiplexing/references/hg38:/hg38 Demuxafy.sif souporcell_pipeline.py \
    -i $BAM \
    -b $BARCODES \
    -f $FASTA \
    -t $THREADS \
    -o $OUT_DIR \
    -k $N \
    --known_genotypes $VCF
       
#Souporcell Summary
singularity exec --bind /omics/odcf/analysis/OE0290_projects/heroes-aya_pools/AG_Thongjuea/Result/10x_multiomics/Demultiplexing_results/demuxafy/RNA/$tool:/$tool Demuxafy.sif bash ${tool}_summary.sh $OUT_DIR/clusters.tsv > /omics/odcf/analysis/OE0290_projects/heroes-aya_pools/AG_Thongjuea/Result/10x_multiomics/Demultiplexing_results/demuxafy/RNA/$tool/$pool/${tool}_summary.tsv

#Correlating Cluster to Donor Reference SNP Genotypes
singularity exec --bind /home/i439h/projects/heroes-aya-pools/AG_Thongjuea/Result/10x_multiomics/Demultiplexing_results/vcf_files:/vcf_files \
        --bind /omics/odcf/analysis/OE0290_projects/heroes-aya_pools/AG_Thongjuea/Result/10x_multiomics/Demultiplexing_results/demuxafy/RNA/$tool:/$tool Demuxafy.sif Assign_Indiv_by_Geno.R \
        -r $VCF -c $OUT_DIR/cluster_genotypes.vcf -o $OUT_DIR

echo "Pool: $pool"

Rscript souporcell_assign.R "$pool"