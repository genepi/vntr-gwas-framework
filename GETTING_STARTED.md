# Getting Started with 1000 Genomes Test Data

UK Biobank data cannot be shared, so this guide runs the pipeline on ten publicly available samples from the [1000 Genomes Project](https://www.internationalgenome.org/) phase 3, using the [30x high-coverage data](https://www.internationalgenome.org/data-portal/data-collection/30x-grch38) (GRCh38).

The guide covers the data preparation part of the pipeline (Steps 1–4). It produces the two inputs needed for GWAS and fine-mapping: a merged VCF with repetitive and non-repetitive *LPA* variants, and per-sample KIV-2 copy numbers.

All steps are reproducible: Steps 2 and 4 run as fixed releases of Nextflow pipelines, which run all tools in Docker containers, and Steps 1 and 3 use tools installed from a single conda environment.

| Step | Description | Runs on test data | Reproducibility |
|---|---|---|---|
| 1 | Extract *LPA*-region reads | ✓ (BAMs provided) | conda |
| 2 | Call KIV-2 VNTR variation | ✓ (output provided) | Nextflow |
| 3 | Combine non-repetitive with repetitive region | ✓ | conda |
| 4 | Estimate KIV-2 copy number | ✓ | Nextflow |
| 5–7 | GWAS, fine-mapping, dosage extraction | UK Biobank only | |

## Setup

All required software is listed in [`scripts/getting_started/environment.yml`](scripts/getting_started/environment.yml). From the repository root, create the environment and a working folder in which all steps are executed:
```
conda env create --file scripts/getting_started/environment.yml
conda activate vntr-getting-started
mkdir getting-started && cd getting-started
```

## Step 1 - Extract *LPA*-region reads

The *LPA*-region BAMs for all ten samples are already provided in [`input/bams/`](input/bams/). They were extracted from the 30x CRAMs as follows (example for HG00265; the CRAM location of each sample is listed in the [sequence index](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/1000G_2504_high_coverage.sequence.index)):
```
samtools view -b -o ../input/bams/HG00265.final.cram.LPA.bam \
  -T http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa \
  http://ftp.sra.ebi.ac.uk/vol1/run/ERR324/ERR3240222/HG00265.final.cram \
  chr6:160530485-160665259
```

## Step 2 - Call KIV-2 VNTR variation

KIV-2 variants are called with [vntr-calling-nf](https://github.com/genepi/vntr-calling-nf) using the signature-sequence approach. We run it on the three European samples (HG00265, HG01685, NA20772). [Docker](https://docs.docker.com/get-docker/) must be installed and running.

Download the single-repeat KIV-2 reference:
```
wget https://raw.githubusercontent.com/genepi/vntr-calling-nf/v0.4.10/reference-data/kiv2.fasta
wget https://raw.githubusercontent.com/genepi/vntr-calling-nf/v0.4.10/reference-data/kiv2.fasta.fai
```

Create `step2.config`:
```
params.project="1000g_eur"
params.input="../input/bams/{HG00265,HG01685,NA20772}.final.cram.LPA.bam"
params.reference="kiv2.fasta"
params.contig="KIV2_6"
params.build="hg38"
params.publish_realigned=true
```

Run the pipeline:
```
nextflow run genepi/vntr-calling-nf -r v0.4.10 -c step2.config -profile docker
```

The VNTR calls are written to `output/1000g_eur/variant_calling/` and the realigned BAMs needed for Step 4 to `output/1000g_eur/realign_fastq/`. Our output of this step is provided with the same structure in [`input/step2-output/`](input/step2-output/), so you can compare your results with it.

## Step 3 - Combine non-repetitive with repetitive region

The VNTR calls are merged with variant calls for the non-repetitive *LPA* region. For 1000 Genomes, these come from the [high-coverage phased panel](https://www.internationalgenome.org/data-portal/data-collection/30x-grch38), which was called from the same 30x CRAMs used in Step 1. The *LPA*-region VCF for all ten samples is already provided in [`input/vcf/`](input/vcf/). It was extracted as follows:
```
bcftools view -r chr6:160530485-160665259 -c 1 \
  -s HG00265,HG00766,HG01685,HG02697,HG03391,HG03673,HG04186,NA18992,NA19087,NA20772 \
  -Oz -o ../input/vcf/lpa_1000g.vcf.gz \
  http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr6.filtered.SNV_INDEL_SV_phased_panel.vcf.gz
```

### 3.1 Convert VNTR results to a VCF file
```
unzip ../scripts/step3/mutserve.zip

# Use sample IDs as names and keep PASS variants only
gzip -dc output/1000g_eur/variant_calling/1000g_eur.txt.gz \
  | sed -E 's/\.[.A-Za-z0-9]*realigned\.bam//g' \
  | awk -F'\t' 'NR==1 || $2=="PASS"' > vntr_filtered.txt

# Rename the KIV-2 contig to 6
sed 's/KIV2_6/6/' kiv2.fasta > kiv2_chr6.fasta

java -jar mutserve.jar create-vcf \
    --input vntr_filtered.txt \
    --output vntr_filtered.vcf.gz \
    --reference kiv2_chr6.fasta
```

### 3.2 Prepare the non-repetitive region
Rename `chr6` to `6` and add an empty DS field, which is filled from the genotypes in the next step:
```
echo "chr6 6" > chr_names.txt
echo '##FORMAT=<ID=DS,Number=1,Type=Float,Description="Genotype dosage">' > ds.hdr

bcftools annotate --rename-chrs chr_names.txt -h ds.hdr -Oz -o region_chr6.vcf.gz ../input/vcf/lpa_1000g.vcf.gz
```

### 3.3 Fix dosages and merge both regions
<!-- TODO (Silvia): provide hg38 coordinates for the KIV-2 exon projection in merge_vntr_nonrep.sh (current values appear to be hg19) and update the text below. -->
The KIV-2 variants are called on a single-repeat reference, where their positions have no meaning on chromosome 6. `merge_vntr_nonrep.sh` projects the variants of both KIV-2 exons onto chromosome 6 coordinates, taking into account that *LPA* lies on the reverse strand. This places repetitive and non-repetitive variants in one VCF that standard GWAS and fine-mapping tools can use directly.
```
sh ../scripts/step3/gt_to_dosage.sh
sh ../scripts/step3/merge_vntr_nonrep.sh vntr_filtered.vcf.gz
sh ../scripts/step3/finalize_dosage.sh
```

The merged VCF is written to `ukb_combined_final_sorted_with_DS_noGT.vcf.gz`. It contains the three European samples.

## Step 4 - Estimate KIV-2 copy number

KIV-2 copy number is estimated from coverage: mean coverage of the KIV-2 exons in the realigned BAMs (Step 2) divided by the coverage of unique *LPA* exons in the original BAMs (Step 1).

Prepare the folders expected by the script and run it:
```
cp -r ../input/bams CRAMS
cp -r output/1000g_eur/realign_fastq realigned
cp -r ../scripts/step4/input input
bash ../scripts/step4/calc_estimates.sh
```

Calculate copy numbers in R (same formula as in [`phenotype.Rmd`](scripts/step4/phenotype.Rmd)):
```r
library(dplyr)
library(tidyr)

read.table("coverage_summary_ukb.txt", header = TRUE, sep = "\t") %>%
  mutate(sample = sub("\\..*", "", BAM)) %>%
  pivot_wider(id_cols = sample, names_from = BED, values_from = SUM) %>%
  mutate(
    cne_kiv2_1 = `kiv2-1.bed` / (1/8 * `exons1.bed`),
    cne_kiv2_2 = `kiv2-2.bed` / (1/8 * `exons2.bed`),
    cne_kiv2 = 1/2 * (cne_kiv2_1 + cne_kiv2_2)
  ) %>%
  filter(!is.na(cne_kiv2))
```

This returns the KIV-2 copy number (`cne_kiv2`) for the three European samples.

## Steps 5–7 - Association and fine-mapping

The remaining steps need Lp(a) measurements and a large cohort, so they cannot be run on the 1000 Genomes test data. In the study, they were run on UK Biobank as described in the [README](README.md):

| Step | Description | Input from this guide |
|---|---|---|
| [5 - GWAS](README.md#step-5---run-combined-gwas-for-lpa-trait) | Genome-wide association for Lp(a) with regenie via [nf-gwas](https://github.com/genepi/nf-gwas), with KIV-2 copy number as covariate | Merged VCF (Step 3), copy numbers (Step 4) |
| [6 - Fine-mapping](README.md#step-6---fine-map-association-signals-using-susie) | SuSiE fine-mapping of the *LPA* locus using an LD matrix from the merged VCF | Merged VCF (Step 3), GWAS results (Step 5) |
| [7 - Dosages](README.md#step-7---extract-dosages-for-credible-set-variants) | Per-sample dosages of credible-set variants for downstream analyses | Fine-mapping results (Step 6) |

## Explore the results

The results of Steps 5 and 6 for UK Biobank are included in this repository:

- **GWAS summary statistics and LD matrices** for European, African and Asian ancestries in [`results/`](results/)
- **[LD Explorer](https://genepi.shinyapps.io/ld-explorer/)**: interactive Manhattan and LD plots of the *LPA* locus for each ancestry
