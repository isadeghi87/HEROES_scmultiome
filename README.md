# HEROES_scmultiome
<<<<<<< HEAD
Codes for analysis of 10x scMultiome including scRNA and scATAC from sarcoma tumors
=======

**Codes for analysis of 10x Genomics Single Cell Multiome (scRNA + scATAC) data from sarcoma tumors**

## Overview

This repository provides a complete computational workflow for processing and analyzing 10x Genomics Single Cell Multiome datasets, integrating both scRNA-seq and scATAC-seq data generated from sarcoma tumor samples. The pipeline covers:

- **Demultiplexing** of pooled samples
- **Alignment & quantification** using Cell Ranger ARC
- **BAM file merging** for pooled lanes or runs
- **Downstream analyses** for scRNA-seq (gene expression) and scATAC-seq (chromatin accessibility)
- **Multi-modal integration** and visualization in Jupyter notebooks

## Features

- **Automated demultiplexing** of pooled 10x scMultiome runs
- **Wrapper scripts** for Cell Ranger ARC to standardize reference usage and parameters
- **Seamless merging** of BAM files for aggregated analyses
- **R-based downstream workflows** using Seurat & Signac for clustering, differential analysis, and gene activity estimation
- **Notebook-driven integration** for exploring correlations between expression and accessibility

## Requirements

The following software and libraries are required:

- **Cell Ranger ARC** (10x Genomics) for initial alignment and quantification
- **Conda** (miniconda or Anaconda)
- **Python** ≥ 3.8
  - `pandas`, `numpy`, `scikit-learn`
- **R** ≥ 4.1
  - `Seurat` (v4), `Signac`, `tidyverse`, `patchwork`, `GenomicRanges`


## Installation

1. **Clone the repository**
   ```bash
   git clone https://github.com/isadeghi87/HEROES_scmultiome.git
   cd HEROES_scmultiome
   ```

2. **Set up Conda environment**
   ```bash
   conda create -n heroes_multiome python=3.9 r-base=4.1
   conda activate heroes_multiome
   # Python deps
   pip install pandas numpy scikit-learn jupyter
   # R deps
   Rscript -e "install.packages(c('Seurat', 'Signac', 'tidyverse', 'patchwork', 'GenomicRanges'))"
   ```

3. **Install Cell Ranger ARC**
   Follow the official guide: https://support.10xgenomics.com/single-cell-multiome/software/downloads

## Directory Structure

```
├── Demultiplexing/             # Scripts to demultiplex pooled scMultiome FASTQs
├── cellRanger_arc/             # Configuration and wrapper for Cell Ranger ARC
├── cellranger_demux_pipeline/  # Combined demultiplex and Cell Ranger ARC pipeline
├── demultiplex_pipeline/       # Alternative demultiplex pipeline (e.g., using bcl2fastq)
├── merge_bams/                 # Scripts to merge BAM files across lanes or samples
├── downstream_analyses/        # R scripts for scRNA and scATAC analyses (Seurat & Signac)
├── integration/                # Jupyter notebooks for multi-modal integration
│   └── .ipynb_checkpoints/     # Auto-saved notebook states
└── README.md                   # This README file
```
>>>>>>> c40d4e2e789bc243d375263f614f667a34e3071a
