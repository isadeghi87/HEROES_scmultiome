#!/bin/bash
source ~/.bash_profile
module load bcftools/1.9
module load vcftools/0.1.16

# here we run vireoSNP tool to demultiplex 10x samples with mitochondrial variants
# Get the pool and ID from command-line arguments
sample_number=$1
ID=$2

sample="pool${sample_number}_$ID"
pool="pool${sample_number}"
OUT_DIR=/home/i439h/projects/heroes-aya-pools/AG_Thongjuea/Result/10x_multiomics/Demultiplexing_results/vireo/results/multiomix/RNA/mito/${sample}
CELL=/home/i439h/projects/heroes-aya-pools/AG_Thongjuea/Result/10x_multiomics/Demultiplexing_results/cellsnp/results/multiomix/RNA/mito/${sample}/cellSNP.cells.vcf.gz
CELL_DATA=/home/i439h/projects/heroes-aya-pools/AG_Thongjuea/Result/10x_multiomics/Demultiplexing_results/cellsnp/results/multiomix/RNA/mito/${sample}/
GT=/home/i439h/projects/heroes-aya-pools/AG_Thongjuea/Result/bulkWGS/mitochondria/vcfs/$pool.merged.vcf.gz

n_donor=4

## create results subfolder for the sample 
if [ ! -d "$OUT_DIR" ]; then
  mkdir -p "$OUT_DIR"
fi

cd $OUT_DIR

vireo -c $CELL_DATA -o $OUT_DIR -d $GT -N $n_donor --genoTag=GT --randSee=12 -p 20 