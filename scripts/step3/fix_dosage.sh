#!/bin/bash

file_no_ds="region_chr6.vcf.gz"
final_file="region_chr6_fixed.vcf.gz"

# Extract variants
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' $file_no_ds > variants.tsv

# Extract DS and GT columns properly
bcftools query -f '[%DS\t]' $file_no_ds | sed 's/\t$//' > ds_only.tsv
bcftools query -f '[%GT\t]' $file_no_ds | sed 's/\t$//' > gt_only.tsv

# Combine
paste variants.tsv ds_only.tsv gt_only.tsv > combined.tsv

# Merge DS and GT
awk '{
    printf "%s\t%s\t%s\t%s", $1,$2,$3,$4;
    n = (NF-4)/2;
    for(i=1;i<=n;i++){
        ds = $(4+i);
        gt = $(4+n+i);
        if(ds=="." || ds=="") {
            # apply GT mapping if ds is missing
            if(gt=="0/0") val=0.0;
            else if(gt=="0/1" || gt=="1/0") val=1.0;
            else if(gt=="1/1") val=2.0;
            else val=".";   # fallback if unexpected
            printf "\t%s", val
        } else {
            printf "\t%s", ds
        }
    }
    printf "\n"
}' combined.tsv > final_unsorted.tsv

rm combined.tsv

sort -k1,1 -k2,2n final_unsorted.tsv > final.tsv
bgzip -c final.tsv > final.tsv.gz
tabix -s 1 -b 2 -e 2 final.tsv.gz

# Create DS header
echo '##FORMAT=<ID=DS,Number=1,Type=Float,Description="DS: overwritten from DS or GT mapping">' > ds_header.txt

# Annotate original VCF with new DS field
bcftools annotate -a final.tsv.gz -h ds_header.txt -c CHROM,POS,REF,ALT,FORMAT/DS -Oz -o $final_file $file_no_ds
