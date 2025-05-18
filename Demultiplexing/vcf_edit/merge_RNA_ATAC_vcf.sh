#!/bin/bash
module load bcftools/1.9
module load htslib/1.8

## VCF files
rna="/home/i439h/projects/pools/AG_Thongjuea/Result/10x_multiomics/Demultiplexing_results/demuxafy/RNA/no_genotype/vireo/pool8/GT_donors.vireo.vcf.gz"
atac="/home/i439h/projects/pools/AG_Thongjuea/Result/10x_multiomics/Demultiplexing_results/demuxafy/ATAC/no_genotype/vireo/pool8/GT_donors.vireo.vcf.gz"
output_vcf="pool8_RNA_ATAC.vcf.gz"

## Working directory
wd="/home/i439h/projects/pools/AG_Thongjuea/Dataset/10x_multiomics/vcfs"
cd $wd
mkdir -p rna_atac_vcfs
cd rna_atac_vcfs

# Create the header file
cat > header.txt << EOL
##fileformat=VCFv4.2
##INFO=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth">
##INFO=<ID=OTH,Number=1,Type=Integer,Description="Other information">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  donor0  donor1  donor2  donor3
EOL

# Decompress VCF files
gunzip -c $rna > rna.vcf
gunzip -c $atac > atac.vcf

# Remove existing headers from VCF files (if any)
grep -v "^##" rna.vcf | grep -v "^#CHROM" > rna_no_header.vcf
grep -v "^##" atac.vcf | grep -v "^#CHROM" > atac_no_header.vcf

# Add new headers to VCF files
cat header.txt rna_no_header.vcf > rna_with_header.vcf
cat header.txt atac_no_header.vcf > atac_with_header.vcf

# Compress the VCF files
bgzip -c rna_with_header.vcf > rna_with_header.vcf.gz
bgzip -c atac_with_header.vcf > atac_with_header.vcf.gz

# Validate the header lines
bcftools view -h rna_with_header.vcf.gz
bcftools view -h atac_with_header.vcf.gz

# Sort the VCF files
bcftools sort -Oz -o rna_with_header_sorted.vcf.gz rna_with_header.vcf.gz
bcftools sort -Oz -o atac_with_header_sorted.vcf.gz atac_with_header.vcf.gz

# Index the VCF files
bcftools index rna_with_header_sorted.vcf.gz
bcftools index atac_with_header_sorted.vcf.gz

# Merge the VCF files
bcftools merge -Oz -o $output_vcf rna_with_header_sorted.vcf.gz atac_with_header_sorted.vcf.gz

# Index the merged VCF file
bcftools index $output_vcf

echo "Merged VCF files saved to $output_vcf"
