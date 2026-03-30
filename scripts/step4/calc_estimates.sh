#!/bin/bash

# Input BED files
BED_FILES=("input/exons1.bed" "input/exons2.bed")
BED_FILES_KIV2=("input/kiv2-1.bed" "input/kiv2-2.bed")
# Directory containing BAM files
BAM_DIR="CRAMS"
BAM_DIR_KIV2="realigned/"
OUTPUT_FILE="coverage_summary_ukb.txt"

# Clear previous output and add header
echo -e "BAM\tBED\tSUM" > "$OUTPUT_FILE"

find "$BAM_DIR" -maxdepth 1 -name "*.bam" -print0 |
while IFS= read -r -d '' BAM; do
    for BED in "${BED_FILES[@]}"; do
        SUM=$(bedtools coverage -a "$BED" -b "$BAM" -mean \
              | awk '{sum += $5} END {print sum+0}')
        echo -e "$(basename "$BAM")\t$(basename "$BED")\t$SUM" >> "$OUTPUT_FILE"
    done
done


# Loop over all BAM files and BED files
find "$BAM_DIR_KIV2" -maxdepth 1 -name "*.bam" -print0 |
while IFS= read -r -d '' BAM; do
    for BED in "${BED_FILES_KIV2[@]}"; do
        SUM=$(bedtools coverage -a "$BED" -b "$BAM" -mean \
              | awk '{sum += $5} END {print sum+0}')
        echo -e "$(basename "$BAM")\t$(basename "$BED")\t$SUM" >> "$OUTPUT_FILE"
    done
done
