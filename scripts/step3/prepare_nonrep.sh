#!/bin/bash
set -e
# Prepare the non-repetitive LPA region: rename chr6 to 6 and add an empty DS field,
# which is filled from the genotypes by gt_to_dosage.sh
# Usage: prepare_nonrep.sh <lpa_region.vcf.gz>
input="$1"

echo "chr6 6" > chr_names.txt
echo '##FORMAT=<ID=DS,Number=1,Type=Float,Description="Genotype dosage">' > ds.hdr

bcftools annotate --rename-chrs chr_names.txt -h ds.hdr -Oz -o region_chr6.vcf.gz "$input"
