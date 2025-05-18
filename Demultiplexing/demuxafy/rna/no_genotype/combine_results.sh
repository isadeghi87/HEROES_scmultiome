
#!/bin/bash

module load anaconda3/2021.05
module load singularity/3.10.5
module load bcftools/1.9

##files
pool="${1:?Pool not specified}"
OUT_DIR=/demuxafy/RNA/no_genotype/combined/$pool
#DEMUXALOT_OUTDIR=/demuxafy/RNA/no_genotype/demuxalot/pool3
SOUPORCELL_OUTDIR=/demuxafy/RNA/no_genotype/souporcell/$pool
VIREO_OUTDIR=/demuxafy/RNA/no_genotype/vireo/$pool
PROJECT_DIR="/home/i439h/projects/pools/AG_Thongjuea"
DEMUX_DIR="${PROJECT_DIR}/Software/10x_multiomics/demuxafy"
cd "$DEMUX_DIR" 

mkdir -p /omics/odcf/analysis/OE0290_projects/heroes-aya_pools/AG_Thongjuea/Result/10x_multiomics/Demultiplexing_results/demuxafy/RNA/no_genotype/combined/$pool

singularity exec --bind /omics/odcf/analysis/OE0290_projects/heroes-aya_pools/AG_Thongjuea/Result/10x_multiomics/Demultiplexing_results/demuxafy:/demuxafy Demuxafy.sif Combine_Results.R \
  -o $OUT_DIR/combined_results.tsv \
  --souporcell $SOUPORCELL_OUTDIR \
  --vireo $VIREO_OUTDIR \
  --method "MajoritySinglet" ## there are other methods that can also be used, please see the help message above for the other options