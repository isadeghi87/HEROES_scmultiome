 #bsub -R "rusage[mem=200G]" -q verylong -n 30 -J 'merge_bam' < merge_bam_v2.sh
 bsub -R "rusage[mem=150G]" -q verylong -n 30 -J 'extratc_reads' < extract_reads.sh