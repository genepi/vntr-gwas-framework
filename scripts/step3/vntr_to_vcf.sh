#!/bin/bash
set -e
# Convert the VNTR calls of vntr-calling-nf (Step 2.1) to a VCF file
# Usage: vntr_to_vcf.sh <vntr_calls.txt.gz> <kiv2.fasta> <output.vcf.gz>
vntr_calls="$1"
reference="$2"
output="$3"
script_dir="$(cd "$(dirname "$0")" && pwd)"

unzip -o -q "$script_dir/mutserve.zip"

# Use sample IDs as names (strip the BAM file suffix) and keep PASS variants only
gzip -dc "$vntr_calls" \
  | sed -E 's/\.[.A-Za-z0-9]*realigned\.bam//g' \
  | awk -F'\t' 'NR==1 || $2=="PASS"' > vntr_filtered.txt

# Rename the KIV-2 contig to 6 (positions are projected onto chr6 in merge_vntr_nonrep.sh)
sed 's/KIV2_6/6/' "$reference" > kiv2_chr6.fasta

"$JAVA_HOME/bin/java" -jar mutserve.jar create-vcf \
    --input vntr_filtered.txt \
    --output "$output" \
    --reference kiv2_chr6.fasta
