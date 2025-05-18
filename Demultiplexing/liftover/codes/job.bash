bsub -R "rusage[mem=150G]" -n 20 -q long < picard_liftover.sh
#bsub -R "rusage[mem=100G]" -n 20 < liftover.sh
#bsub -R "rusage[mem=100G]" -n 20 < mergeVcf.sh