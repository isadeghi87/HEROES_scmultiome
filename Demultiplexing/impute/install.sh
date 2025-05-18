#!/bin/bash
module load anaconda3/2021.05
wd=/home/i439h/projects/pools/AG_Thongjuea/Software/10x_multiomics/impute

# Create a Conda environment in a specific directory
conda create --prefix $wd/impute_env python=3.9
source activate $wd/impute_env 
conda install -c bioconda impute2

