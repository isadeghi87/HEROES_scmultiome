#!/bin/bash
bsub -R "rusage[mem=100G]" -q long -n 10 -J merge <merge_bulkWGS_pools.sh
#bsub -R "rusage[mem=100G]" -q long -n 2 -J filter <filter_qual3000_dp100.sh 
