set -e
echo "get lpa non-rep region"
#bcftools filter /mnt/genepi-biobank/data/gwas/ukbb/imputed/vcfs/ukb_imp_chr6_v3.vcf.gz  --regions 6:160952514-161033863,6:161038409-161085307 -Ov -o ukb_lpa_nonrep.vcf.gz
# THIS FILE HAS BEEN CREATED BY SETTING THE GT value to DS for GENOTYPED ONLY VARIANTS 
cp region_chr6_fixed.vcf.gz ukb_lpa_nonrep.vcf.gz
tabix ukb_lpa_nonrep.vcf.gz
echo "extract exome regions"
bcftools norm -m -any ukb_rap_renamed_filtered.vcf.gz  -o 2022-10-21-ukbb-kiv2-annoated.renamed.passed.norm.vcf.gz -Oz
tabix -p vcf -f 2022-10-21-ukbb-kiv2-annoated.renamed.passed.norm.vcf.gz
bcftools filter 2022-10-21-ukbb-kiv2-annoated.renamed.passed.norm.vcf.gz --regions 6:481-840 -Ov -o ukb_kiv2.6_exon1.vcf.gz
bcftools filter 2022-10-21-ukbb-kiv2-annoated.renamed.passed.norm.vcf.gz --regions 6:4644-5025 -Ov -o ukb_kiv2.6_exon2.vcf.gz

echo "set rsID with chrom + pos" 
bcftools annotate --set-id +'%CHROM\:%POS\:%REF\:%ALT' ukb_kiv2.6_exon1.vcf.gz -Oz -o ukb_kiv2.6_exon1.ids.vcf.gz
bcftools annotate --set-id +'%CHROM\:%POS\:%REF\:%ALT' ukb_kiv2.6_exon2.vcf.gz -Oz -o ukb_kiv2.6_exon2.ids.vcf.gz

bcftools sort -Oz -o ukb_kiv2.6_exon1.ids.sorted.vcf.gz ukb_kiv2.6_exon1.ids.vcf.gz
bcftools sort -Oz -o ukb_kiv2.6_exon2.ids.sorted.vcf.gz ukb_kiv2.6_exon2.ids.vcf.gz
tabix -p vcf ukb_kiv2.6_exon1.ids.sorted.vcf.gz
tabix -p vcf ukb_kiv2.6_exon2.ids.sorted.vcf.gz

bcftools concat -a -Oz -o ukb_kiv2.6_combined.vcf.gz \
    ukb_kiv2.6_exon1.ids.sorted.vcf.gz \
    ukb_kiv2.6_exon2.ids.sorted.vcf.gz
tabix -p vcf -f ukb_kiv2.6_combined.vcf.gz

zcat ukb_kiv2.6_combined.vcf.gz | \
awk -F'\t' 'BEGIN{OFS="\t"}
/^#/ {print; next}
{
    split($9, fmt, ":")           # FORMAT keys
    for(s=10; s<=NF; s++){        # for each sample column
        split($s, vals, ":")
        for(i=1;i<=length(fmt);i++) f[fmt[i]]=i
        if(vals[f["GT"]]=="0/0"){
            for(i=1;i<=length(vals);i++){
                if(i!=f["GT"]) vals[i]="."    # set all except GT to "."
            }
        }
        $s=vals[1]
        for(i=2;i<=length(vals);i++) $s=$s ":" vals[i]
    }
    print
}' > ukb_kiv2.6_combined_fixed.vcf

bgzip ukb_kiv2.6_combined_fixed.vcf
tabix -p vcf -f ukb_kiv2.6_combined_fixed.vcf.gz

zcat ukb_kiv2.6_combined_fixed.vcf.gz | \
awk 'BEGIN {
    OFS = "\t"
    # exon1 mapping
    start1 = 481
    end1   = 840
    map_start1 = 161038408
    map_end1   = 161038049

    # exon2 mapping
    start2 = 4644
    end2   = 5025
    map_start2 = 161034245
    map_end2   = 161033864
}
# keep headers unchanged
/^#/ { print; next }

{
    # exon1 transform
    if ($2 >= start1 && $2 <= end1) {
        scale1 = (map_start1 - map_end1) / (end1 - start1)
        new_pos = map_start1 - (($2 - start1) * scale1)
        $2 = int(new_pos)
    }

    # exon2 transform
    else if ($2 >= start2 && $2 <= end2) {
        scale2 = (map_start2 - map_end2) / (end2 - start2)
        new_pos = map_start2 - (($2 - start2) * scale2)
        $2 = int(new_pos)
    }

    print
}' > ukb_kiv2.6_combined_remapped.vcf

bgzip ukb_kiv2.6_combined_remapped.vcf

bcftools sort -Oz -o ukb_kiv2.6_combined_remapped_sorted.vcf.gz ukb_kiv2.6_combined_remapped.vcf.gz
tabix -p vcf -f ukb_kiv2.6_combined_remapped_sorted.vcf.gz

echo "filter by samples"
bcftools query -l ukb_kiv2.6_combined_remapped_sorted.vcf.gz > samples_repetitive.txt
bcftools query -l ukb_lpa_nonrep.vcf.gz > samples_nonrepetitive.txt
grep -Fxf samples_nonrepetitive.txt samples_repetitive.txt > common_ids.txt
bcftools view -S common_ids.txt -Oz -o ukb_kiv2.6_combined_remapped_sorted_filtered.vcf.gz ukb_kiv2.6_combined_remapped_sorted.vcf.gz
bcftools view -S common_ids.txt -Oz -o ukb_lpa_nonrep_filtered.vcf.gz ukb_lpa_nonrep.vcf.gz
tabix -p vcf -f ukb_kiv2.6_combined_remapped_sorted_filtered.vcf.gz
tabix -p vcf -f ukb_lpa_nonrep_filtered.vcf.gz

echo "concat regions"
bcftools concat -a -Oz -o ukb_combined_final.vcf.gz \
    ukb_kiv2.6_combined_remapped_sorted_filtered.vcf.gz \
    ukb_lpa_nonrep_filtered.vcf.gz

bcftools sort ukb_combined_final.vcf.gz -Oz -o ukb_combined_final_sorted.vcf.gz
tabix -p vcf -f ukb_combined_final_sorted.vcf.gz
