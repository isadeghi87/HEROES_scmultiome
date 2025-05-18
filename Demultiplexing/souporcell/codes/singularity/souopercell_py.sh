module load anaconda3/2021.05
source activate /home/i439h/.conda/envs/souporcell
 source ~/.bash_profile

sample='pool1_1016711'
BAM=/omics/odcf/analysis/OE0536_projects/heroes-aya/scRNA-Seq/Demultiplexing_results/souporcell/input_files/${sample}/${sample}_genome_bam.bam
BARCODE=/omics/odcf/analysis/OE0536_projects/heroes-aya/scRNA-Seq/Demultiplexing_results/souporcell/input_files/${sample}/barcodes.tsv
FASTA=/omics/odcf/analysis/OE0536_projects/heroes-aya/scRNA-Seq/CellRanger_v7.1_results/References/10x/refdata-gex-GRCh38-2020-A/fasta/genome.fa
OUT_DIR=/omics/odcf/analysis/OE0536_projects/heroes-aya/scRNA-Seq/Demultiplexing_results/souporcell/results/${sample}/

souporcell_pipeline.py -i $BAM \
    -b $BARCODE \
    -f $FASTA \
    -t 22 \
    -o $OUT_DIR \
    -k 4