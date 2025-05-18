bsub -R "rusage[mem=100G]" -q long -n 10 < souporcell_pool1.sh
bsub -R "rusage[mem=100G]" -q long -n 10 < souporcell_pool2.sh
bsub -R "rusage[mem=100G]" -q long -n 10 < souporcell_pool3.sh
bsub -R "rusage[mem=100G]" -q long -n 10 < souporcell_pool4.sh


bsub -R "rusage[mem=200G]" -q long -n 10 < pool3_multiomix_RNA_WGS_BJ.sh
bsub -R "rusage[mem=200G]" -q long -n 10 < pool4_multiomix_RNA_WGS_BJ.sh