#!/bin/bash
## liftover using python module wrapper for ucsc liftover https://github.com/single-cell-genetics/cellSNP/tree/master/liftOver

chain=/home/i439h/projects/heroes-aya/scRNA-Seq/Demultiplexing_results/references/liftover/hg19ToHg38.over.chain 

#python liftOver_vcf.py -c $chain -i my_input.hg19.vcf.gz -o my_output.hg38.vcf.gz
python liftOver_vcf.py -c $chain -i my_input.hg19.vcf.gz -o my_output.hg38.vcf.gz