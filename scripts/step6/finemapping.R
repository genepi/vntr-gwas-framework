## -----------------------------
## Who / When / What
## -----------------------------
# Silvia Di Maio (based on Johanna F. Schachtl-Riess)
# 19.12.2025
# Remove covariate effects from genetic data (X) for VNTR AFR project

## -----------------------------
## Packages (minimal)
## -----------------------------
library(Matrix)      # matrix operations, crossprod, nearPD
library(data.table)  # fread, fwrite
library(vcfR)        # read VCF, extract DS
library(tidyverse)   # dplyr, tibble, ggplot2
library(stats)       # lm, residuals, cor
library(graphics)    # plot, abline
library(susieR)      # SuSiE fine-mapping
library(R.utils)     # utility functions
library(dplyr)

## -----------------------------

## -----------------------------
## Load VCF and fix duplicate IDs
## -----------------------------
cat("Reading VCF...\n")
vcf <- read.vcfR("ukb_kiv2_estimates_final_sorted_with_DS_noGT_afr_filtered.vcf.gz")

cat("Creating unique IDs using CHROM_POS_REF_ALT...\n")
safe_alt <- gsub(",", "_", vcf@fix[, "ALT"])  # handle multi-allelic ALT
unique_ids <- paste(vcf@fix[, "CHROM"], vcf@fix[, "POS"], vcf@fix[, "REF"], safe_alt, sep="_")

# Overwrite ID column
vcf@fix[, "ID"] <- unique_ids
rownames(vcf@fix) <- unique_ids
rownames(vcf@gt) <- unique_ids

# Sanity check
stopifnot(!any(duplicated(vcf@fix[, "ID"])))

cat("Extracting DS matrix...\n")
ds_matrix <- extract.gt(vcf, element = "DS", as.numeric = TRUE)

# Transpose: samples x SNPs
X <- t(ds_matrix)
stopifnot(is.matrix(X), is.numeric(X))
cat("DS matrix dimensions:", dim(X), "\n")
cat("Number of NAs:", sum(is.na(X)), "\n")

## -----------------------------
## Load covariates
## -----------------------------
covariates <- fread("input/phenotype_ukb_estimates_afr_covariates.txt", header = TRUE)
nrow(covariates)
# Remove unnecessary columns
covariates$FID <- as.factor(covariates$FID)
covariates <- covariates %>% dplyr::select(-IID, -lpa_man, -exons1_sum, -exons2_sum, 
                                           -ancestry, -kiv2_1, -kiv2_2, -cne_kiv2_1, -cne_kiv2_2, -inv_norm_BL_ibk_lpa)

# Convert to matrix with FID as rownames
cov_matrix <- covariates %>%
  tibble::column_to_rownames("FID") %>%
  filter(complete.cases(.)) %>%
  as.matrix()
stopifnot(!any(is.na(cov_matrix)))

# check for NAs
any(is.na(cov_matrix))

# check dimension
dim(cov_matrix)


## -----------------------------
## Align samples
## -----------------------------
common_names <- intersect(rownames(X), rownames(cov_matrix))
X <- X[common_names, , drop = FALSE]
cov_matrix <- cov_matrix[common_names, , drop = FALSE]

# check if individuals are in the same order
identical(rownames(cov_matrix), 
          rownames(X))

## -----------------------------
## Remove covariate effects
## -----------------------------
remove.covariate.effects <- function(X, Z) {
  if (any(Z[,1] != 1)) Z <- cbind(1, Z)
  A <- forceSymmetric(crossprod(Z))
  SZX <- solve(A, t(Z) %*% X)
  X <- X - Z %*% SZX
  list(X = X, SZX = SZX)
}

cat("Removing covariate effects...\n")
out <- remove.covariate.effects(X, cov_matrix)

## -----------------------------
## Sanity check
## -----------------------------
X_resid_5 <- out$X[,5]
lm_resid_5 <- residuals(lm(X[,5] ~ cov_matrix))
cat("Correlation of residuals:", cor(X_resid_5, lm_resid_5), "\n")

plot(X_resid_5, lm_resid_5,
     main = "Residual Comparison for X[,5]",
     xlab = "Manual Residuals", ylab = "lm() Residuals")
abline(0, 1, col = "red")

## -----------------------------
## Correlation matrix and export
## -----------------------------
# Ensure base numeric matrix for cor()
X_residual <- as.matrix(out$X)
storage.mode(X_residual) <- "numeric"

cor_R <- cor(X_residual)
cat("Correlation matrix dims:", dim(cor_R), "\n")
cat("Range:", min(cor_R), "to", max(cor_R), "\n")

write.table(cor_R, file = "UKB_afr_ld_residuals.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

## STEP 2Folder setup =======

dir.create("output", showWarnings = FALSE)

## Data ===========

## READ IN DATA AND CALCULATE Z-Scores
ld_matrix <- as.matrix(fread("UKB_afr_ld_residuals.txt", header = TRUE))
dim(ld_matrix)
regenie <- fread("lpa_man.regenie_afr.gz", header = TRUE)
betahat <- regenie$BETA
sebetahat <- regenie$SE
z_scores <- betahat / sebetahat

## LD matrix of residuals after covariate adjustment ===============

### QC =============

# Check if LD matrix is symmetric and PSD
if (!isSymmetric(ld_matrix)) {
  warning("LD matrix is not symmetric. Adjusting with nearPD.")
}
if (any(eigen(ld_matrix)$values < 0)) {
  warning("LD matrix is not positive semidefinite. Adjusting with nearPD.")
}


#ensure that the linkage disequilibrium (LD) matrix is positive semidefinite (PSD)
ld_matrix <- as.matrix(nearPD(ld_matrix)$mat)

# ensure that variants are in the correct order
dim(ld_matrix)
nrow(regenie)

identical(rownames(ld_matrix), 
          regenie$ID)

# Output first 10 row names of ld_matrix
head(rownames(ld_matrix), 10)

# Output first 10 IDs from regenie
head(regenie$GENPOS, 10)

# remove names from matrix
dimnames(ld_matrix) <- NULL

### Diagnostics ===================
# link: https://stephenslab.github.io/susieR/articles/susierss_diagnostic.html
# goal: assessing consistency of the summary statistics and the reference LD matrix

# calculate lambda
# should be low (rule of thumb: below 0.1?)
lambda = estimate_s_rss(z_scores, 
                        ld_matrix, 
                        n=5726)
lambda

# Compute Distribution of z-scores of Variant j Given Other z-scores, and Detect Possible Allele Switch Issue
# should be on one line if from the same data, if from a reference panel should be close around line, extreme outliers can point to allele switches
condz_in = kriging_rss(z_scores, 
                       ld_matrix, 
                       n=5726)
condz_in$plot

# Add lambda annotation to the plot
annotated_plot <- condz_in$plot +
  annotate("text", 
           x = Inf, y = -Inf, 
           label = paste("Lambda:", format(lambda, digits = 4)), 
           hjust = 1.1, vjust = -1.2, 
           size = 4, fontface = "bold")

# Save the plot with the annotation
ggsave(filename = "output/susie_diagnostics_plot_refined_estimates_adjust_afr.png",
       plot = annotated_plot,
       width = 7, height = 5, dpi = 300)


### SUSIE =================
# 

susie_out_resid <- susie_rss(
  z = z_scores,
  R = ld_matrix,
  n = 5726, 
  L = 50,
  estimate_residual_variance = TRUE,
  check_input = TRUE,  # Enable input checks,
  max_iter = 100,
  refine = TRUE,
)

#### Checks ==========

# should converge else increase max_iter
susie_out_resid$converged

#### Results ===============
summary_susie_out_resid <- summary(susie_out_resid)

summary_susie_out_resid$cs

susie_out_resid$sets$purity

##### CREDIBLE SETS ===========

credible_sets <- susie_out_resid$sets$cs
pip <- susie_out_resid$pip

credible_set_tables_with_ld <- lapply(credible_sets, function(indices) {
  
  ld_subset <- ld_matrix[indices, indices]
  
  if (length(indices) == 1) {
    ld_values <- as.character(round(ld_subset, 2))
  } else if (is.null(dim(ld_subset)) || nrow(ld_subset) == 0 || ncol(ld_subset) == 0) {
    warning("LD subset is empty or invalid for credible set.")
    return(NULL)
  } else {
    ld_values <- apply(ld_subset, 1, function(row) paste(round(row, 2), collapse = ","))
  }
  
  regenie %>%
    slice(indices) %>%
    mutate(PIP = pip[indices]) %>%
    mutate(LD = ld_values) %>%
    dplyr::select(CHROM, GENPOS, ALLELE0, ALLELE1, A1FREQ, LOG10P, 
                  BETA, RSID, N, PIP, LD, ID)
})

# Combine all credible sets into one data.frame
combined_df <- rbindlist(
  lapply(seq_along(credible_set_tables_with_ld), function(i) {
    dt <- credible_set_tables_with_ld[[i]]
    dt[, layer := sub("^L", "", names(credible_set_tables_with_ld)[i])]  # add layer info
    dt
  }),
  use.names = TRUE,
  fill = TRUE
)

fwrite(combined_df, file = "output/afr_credible_sets.txt", sep = "\t")

write.table(combined_df %>% dplyr::select(GENPOS),
            "output/afr_credible_sets_pos.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)

