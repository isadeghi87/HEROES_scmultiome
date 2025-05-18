
bsub -R "rusage[mem=200G]" -q long -n 10 -J pool3_1000genome < pool3_multiomix_RNA_1000genome.sh
#bsub -R "rusage[mem=200G]" -q long -n 10 < pool4_multiomix_RNA_WGS_BJ.sh
