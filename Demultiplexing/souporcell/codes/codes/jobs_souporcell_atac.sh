bsub -R "rusage[mem=100G]" -q verylong -n 10 < pool1_multiomix_atac.sh
bsub -R "rusage[mem=100G]" -q verylong -n 10 < pool2_multiomix_atac.sh
bsub -R "rusage[mem=100G]" -q verylong -n 10 < pool3_multiomix_atac.sh
bsub -R "rusage[mem=100G]" -q verylong -n 10 < pool4_multiomix_atac.sh