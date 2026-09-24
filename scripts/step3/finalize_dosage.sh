file_no_ds="ukb_combined_final_sorted.vcf.gz"
file_ds="ukb_combined_final_sorted_with_DS.vcf.gz"
final="ukb_combined_final_sorted_with_DS_noGT.vcf.gz"
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' $file_no_ds > variants.tsv
bcftools query -f '[%DS\t]\n' $file_no_ds > ds_only.tsv
bcftools query -f '[%AF\t]\n' $file_no_ds > af_only.tsv
wc -l variants.tsv
wc -l ds_only.tsv
wc -l af_only.tsv
paste variants.tsv ds_only.tsv af_only.tsv > combined.tsv

cat combined.tsv | awk '{
    printf "%s\t%s\t%s\t%s", $1, $2, $3, $4;
    n = (NF - 4) / 2;
    for (i = 5; i <= 4 + n; i++) {
        ds = $i;
        af = $(i + n);
        if (ds != ".") {
            printf "\t%s", ds;
        } else if (af != ".") {
            printf "\t%.3f", af + 1;
        } else {
            printf "\t0";
        }
    }
    printf "\n";
}'  > final.tsv


bgzip -c final.tsv > final.tsv.gz
tabix -s 1 -b 2 -e 2 final.tsv.gz
echo '##FORMAT=<ID=DS,Number=1,Type=Float,Description="DS: overwritten from DS or AF+1">' > ds_header.txt

bcftools annotate -a final.tsv.gz   -h ds_header.txt   -c CHROM,POS,REF,ALT,FORMAT/DS   -Oz -o $file_ds $file_no_ds
tabix $file_ds

echo "FILTER"
# this is because when we use NORM the AF gets included in BOTH lines. If we dont add that we would se 4733 twice. So we need a better way than this. Here we filter rare variants (multiallelic now biallelic)
# AC 1000 means 0.5%
#bcftools view -i 'AC>=1000' "$file_ds" \
#  | bcftools annotate -x FORMAT/GT -Oz -o "$final"

bcftools annotate -x FORMAT/GT "$file_ds" -Oz -o "$final"
