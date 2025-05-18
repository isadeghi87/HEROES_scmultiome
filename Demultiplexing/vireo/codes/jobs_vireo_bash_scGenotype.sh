bsub -R "rusage[mem=100G]" -q long -n 20 < vireo_pool1_scGenotype.sh
bsub -R "rusage[mem=100G]" -q long -n 20 < vireo_pool2_scGenotype.sh
bsub -R "rusage[mem=100G]" -q long -n 20 < vireo_pool3_scGenotype.sh
bsub -R "rusage[mem=100G]" -q long -n 20 < vireo_pool4_scGenotype.sh