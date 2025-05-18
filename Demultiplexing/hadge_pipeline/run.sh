## running nexflow nf-core/methylseq pipeline for HEREOES-AYA OE0290 whole genome biosulfite sequencing data 
# load nextflow module 
module load anaconda3/2021.05
source activate env_nf
module load singularity/3.10.5
module load openjdk/17u35

##Directories and files
## codes 
codes=/omics/odcf/analysis/OE0290_projects/heroes-aya_pools/AG_Thongjuea/Result/10x_multiomics/Demultiplexing_results/hadge_pipeline/hadge
#export NXF_SINGULARITY_CACHEDIR=/omics/odcf/analysis/OE0290_projects/pediatric_tumor/whole_genome_biosulfite_sequencing/codes/nextflow/work/singularity
#export NXF_OPTS="-Dleveldb.mmap=false"
fasta=/omics/odcf/analysis/OE0290_projects/heroes-aya_pools/temp_analysis/scRNA-Seq/CellRanger_v7.1_results/References/10x/refdata-gex-GRCh38-2020-A/fasta/genome.fa
fasta_index=/omics/odcf/analysis/OE0290_projects/heroes-aya_pools/temp_analysis/scRNA-Seq/CellRanger_v7.1_results/References/10x/refdata-gex-GRCh38-2020-A/fasta/genome.fa.fai
samples=/home/i439h/projects/heroes-aya-pools/AG_Thongjuea/Result/10x_multiomics/Demultiplexing_results/hadge_pipeline/samples.csv
outdir=/home/i439h/projects/heroes-aya-pools/AG_Thongjuea/Result/10x_multiomics/Demultiplexing_results/hadge_results

cd $codes
nextflow main.nf --multi_input $samples --outdir $outdir -profile conda --fasta $fasta --fasta_index $fasta_index \
    --mode genetic --max_cpus 30 --max_memory 200.GB -bg --generate_anndata False \
    --common_variants_cellsnp ${vcf_donor} --match_donor False

