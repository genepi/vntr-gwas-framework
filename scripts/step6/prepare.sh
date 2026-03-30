set - e
pops=("afr")
for pop in "${pops[@]}"; do
filename_vcf="ukb_kiv2_estimates_final_sorted_with_DS_noGT_${pop}.vcf.gz"
filename_regenie="lpa_man.regenie_${pop}.gz"
filename_out="ukb_kiv2_estimates_final_sorted_with_DS_noGT_${pop}_filtered.vcf.gz"

rm $filename_out
rm $filename_regenie
rm $yfilename_vcf
cp input/lpa_man.regenie.gz $filename_regenie
cp input/ukb_combined_final_sorted_with_DS_noGT_afr.vcf.gz  $filename_vcf

zcat $filename_regenie | awk -F'\t' 'BEGIN{OFS="\t"} NR>1 {print $1, $2, $4, $5}' > regenie.snps

echo #REGENIE
wc -l regenie.snps

# Prepare an output VCF with header
bcftools view -h $filename_vcf > ${filename_out%.gz}

while IFS=$'\t' read -r chrom pos ref alt; do
    # Extract the exact variant (multi-allelic-safe)
    bcftools view -H -i "POS==$pos && REF==\"$ref\" && ALT==\"$alt\"" \
        $filename_vcf >> ${filename_out%.gz}
done < regenie.snps

bgzip -c ${filename_out%.gz} > $filename_out
tabix -p vcf $filename_out

echo #VCF
bcftools view -H $filename_out | wc -l

rm ${filename_out%.gz}
rm $filename_vcf
#rm $filename_regenie
done
