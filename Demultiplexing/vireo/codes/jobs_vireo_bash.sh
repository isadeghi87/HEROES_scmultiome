bsub -R "rusage[mem=100G]" -q long -n 20 < vireo_pool1.sh
bsub -R "rusage[mem=100G]" -q long -n 20 < vireo_pool2.sh
bsub -R "rusage[mem=100G]" -q long -n 20 < vireo_pool3.sh
bsub -R "rusage[mem=100G]" -q long -n 20 < vireo_pool4.sh

bsub -R "rusage[mem=200G]" -q long -n 20 -J vireo_pool3 < vireo_pool3_multiomics_WGS_BJ.sh
bsub -R "rusage[mem=200G]" -q long -n 20 -J vireo_pool4 < vireo_pool4_multiomics_WGS_BJ.sh

bsub -R "rusage[mem=200G]" -q long -n 20 -J vireo_pool3_RNA < vireo_pool3_scRNA_WGS_BJ.sh
bsub -R "rusage[mem=200G]" -q long -n 20 -J vireo_pool4_RNA < vireo_pool4_scRNA_WGS_BJ.sh