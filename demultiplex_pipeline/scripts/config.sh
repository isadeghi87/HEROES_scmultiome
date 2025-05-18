#!/bin/bash

# General directories; change directories to match your paths
OUT_DIR="/home/i439h/projects/heroes-aya-pools/AG_Thongjuea/Code/10x_multiomics/Sadeghi/demultiplex_pipeline/output"
#VCF_INPUT_DIR="/path/to/vcf/input"
VCF_OUTPUT_DIR="$OUT_DIR/vcf/output"
CELLSNP_OUTPUT_DIR="$OUT_DIR/cellsnp/output"
VIREO_OUTPUT_DIR="$OUT_DIR/vireo/output"
SOUPORCELL_OUTPUT_DIR="$OUT_DIR/souporcell/output"
SAMPLE_LIST_DIR="$OUT_DIR/sample_lists"


# Path to the samplesheet
SAMPLESHEET_PATH="$OUT_DIR/samplesheet.tsv"

# Tool-specific configuration
BCFTOOLS_MODULE="bcftools/1.9"
HTSLIB_MODULE="htslib/1.8"
ANACONDA_MODULE="anaconda3/2021.05"
CELLSNP_CONDA_ENV="/home/i439h/projects/heroes-aya-pools/temp_analysis/tools/conda/cellsnp"
VCFTOOLS_MODULE="vcftools/0.1.16"

# Number of donors for Vireo
n_donor=4

# Job scheduler settings (adjust based on your HPC environment)
PIPELINE_QUEUE_NAME='verylong' 
QUEUE_NAME='long' 
FILTER_JOB_MEM="100G"  # Adjust as needed for filter_vcf step
FILTER_JOB_CPUS=20

MERGE_JOB_MEM="100G"  # Adjust as needed for merge_vcf step
MERGE_JOB_CPUS=10

CELLSNP_JOB_MEM="100G"  # Adjust as needed for cellsnp step
CELLSNP_JOB_CPUS=10

VIREO_JOB_MEM="100G"  # Adjust as needed for vireo step
VIREO_JOB_CPUS=10

# SOUPORCELL_JOB_MEM="100G"  # Adjust as needed for souporcell step
# SOUPORCELL_JOB_CPUS=10

# Other global settings
