
#!/bin/bash

module load anaconda3/2021.05
module load singularity/3.10.5
module load bcftools/1.9

##files
OUT_DIR=/demuxafy/RNA/combined
DEMUXALOT_OUTDIR=/demuxafy/RNA/demuxalot/pool3
SOUPORCELL_OUTDIR=/demuxafy/RNA/souporcell/pool3
VIREO_OUTDIR=/demuxafy/RNA/vireo/pool3

mkdir -p /omics/odcf/analysis/OE0290_projects/heroes-aya_pools/AG_Thongjuea/Result/10x_multiomics/Demultiplexing_results/demuxafy/RNA/combined

singularity exec --bind /omics/odcf/analysis/OE0290_projects/heroes-aya_pools/AG_Thongjuea/Result/10x_multiomics/Demultiplexing_results/demuxafy:/demuxafy Demuxafy.sif Combine_Results.R \
  -o $OUT_DIR/combined_results.tsv \
  --demuxalot $DEM  UXALOT_OUTDIR \
  --souporcell $SOUPORCELL_OUTDIR \
  --vireo $VIREO_OUTDIR \
  --method "MajoritySinglet" ## there are other methods that can also be used, please see the help message above for the other options