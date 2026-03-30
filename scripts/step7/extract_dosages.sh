#!/bin/bash
set -euo pipefail

VCF="input/ukb_kiv2_estimates_final_sorted_with_DS_noGT_afr_filtered.vcf.gz"
OUT="snps_dosages_estimates_afr.csv"
POS_FILE="input/afr_credible_sets_pos.txt"


# Remove old temp files if any
rm -f *.csv *.vcf.gz

# Step 1: Extract SNPs and convert to CSV
while read -r pos; do
    # Extract SNP from VCF
    bcftools view -r "6:${pos}" "$VCF" -Oz -o "${pos}.vcf.gz"

    # Convert VCF to CSV using genomic-utils
    java -jar genomic-utils.jar vcf-to-csv-transpose \
        --input "${pos}.vcf.gz" \
        --output "${pos}.csv" \
        --genotypes DS
done < "$POS_FILE"

# Step 2: Merge CSVs
csv_files=($(ls *.csv | sort -V))  # sort numerically by SNP position

# Initialize output with first column of first CSV
cut -d, -f1 "${csv_files[0]}" > "$OUT"

# Loop over all CSVs to append second columns
for f in "${csv_files[@]}"; do
    cut -d, -f2 "$f" > temp_col.csv
    paste -d, "$OUT" temp_col.csv > merged_tmp.csv
    mv merged_tmp.csv "$OUT"
done

# Step 3: Clean up
rm -f temp_col.csv
rm -f *.vcf.gz
rm -f 1*.csv

echo "Merged CSV created at $OUT"
