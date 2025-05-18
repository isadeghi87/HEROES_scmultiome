#!/bin/bash
# load modules
module load anaconda3/2021.05
source activate /home/i439h/.conda/envs/souporcell
module load samtools/1.6
module load bcftools/1.9
module load bedtools/2.24.0
module load  htslib/1.4.1
module load freebayes/1.3.4
module load vartrix/1.1.22
module load minimap2/2.24
source ~/.bash_profile

#module load singularity/3.10.5
# here we run souporcell tool to demultiplex 10x scRNAseq samples without genotyping data
sample='pool4_1016717'
pool='Pool4'
OUT_DIR=/omics/odcf/analysis/OE0536_projects/heroes-aya/scRNA-Seq/Demultiplexing_results/souporcell/results/${sample}/
BAM=/omics/odcf/analysis/OE0536_projects/heroes-aya/scRNA-Seq/Demultiplexing_results/souporcell/input_files/${sample}/${sample}_genome_bam.bam
BARCODE=/omics/odcf/analysis/OE0536_projects/heroes-aya/scRNA-Seq/Demultiplexing_results/souporcell/input_files/${sample}/${sample}_barcodes.tsv.gz
FASTA=/omics/odcf/analysis/OE0536_projects/heroes-aya/scRNA-Seq/CellRanger_v7.1_results/References/10x/refdata-gex-GRCh38-2020-A/fasta/genome.fa
GT=/omics/odcf/analysis/OE0536_projects/heroes-aya/scRNA-Seq/Demultiplexing_results/liftover/vcf_files/$pool.vcf.gz

## create results subfolder for the sample 
if [ ! -d "$OUT_DIR" ]; then
  mkdir -p "$OUT_DIR"
fi

cd $OUT_DIR
cp $GT ./
gunzip $pool.vcf.gz
souporcell_pipeline.py \
    -i $BAM \
    -f $FASTA \
    -t 20 \
    -o $OUT_DIR \
    -b $BARCODE \
    --umi_tag UB \
    --cell_tag CB \
    --known_genotypes $pool.vcf \
    --known_genotypes_sample_names I020_004_metastasis011-01_control001-01 I023_047_tumor091-01_control001-01 I010_022_tumor011-01_control001-01 I005_016_metastasis011-01_control001-02 \
    -k 4