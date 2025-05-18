#!/bin/bash

# Check for the correct number of arguments
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 input_bam output_bed"
    exit 1
fi

input_bam=$1
output_bed=$2

# Extract cell barcodes from BAM using samtools and append to BED
samtools view -h $input_bam | awk -v OFS="\t" '
BEGIN { FS="\t" }
{
    if ($1 ~ /^@/) {
        # Skip header lines
        next
    }
    # Extract cell barcode from CB:Z tag
    for (i = 12; i <= NF; i++) {
        if ($i ~ /^CB:Z:/) {
            cell_barcode = substr($i, 6)
            break
        }
    }
    # Output BED line with cell barcode
    print $3, $4, $4 + length($10), $1, 0, $2, cell_barcode
}' | sort -k1,1 -k2,2n > $output_bed

echo "Finished processing $input_bam to $output_bed"
