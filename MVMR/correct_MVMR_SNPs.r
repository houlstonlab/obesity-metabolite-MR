###
# Script that corrects output from mv_extract_exposures_local by
# removing palindromic SNPs that cannot be inferred.
# This has now been fixed in the version 0.6.12 of TwoSampleMR.
#
# args[1] = cancer name
###

args <- commandArgs(trailingOnly = TRUE)

outdir <- "/path/to/your/output/directory/"
correct_harmonised <- read.table(paste0(outdir, "collated_harmonised.tsv"), header = TRUE, row.names = 1)
incorrect_harmonised <- readRDS(paste0(outdir, "cancer_harmonised.rds"))

mask <- rownames(incorrect_harmonised$exposure_beta) %in% rownames(correct_harmonised)

incorrect_harmonised$exposure_beta <- incorrect_harmonised$exposure_beta[mask, ]
incorrect_harmonised$exposure_se <- incorrect_harmonised$exposure_se[mask, ]
incorrect_harmonised$exposure_pval <- incorrect_harmonised$exposure_pval[mask, ]
incorrect_harmonised$outcome_beta <- incorrect_harmonised$outcome_beta[mask]
incorrect_harmonised$outcome_se <- incorrect_harmonised$outcome_se[mask]
incorrect_harmonised$outcome_pval <- incorrect_harmonised$outcome_pval[mask]

saveRDS(incorrect_harmonised, paste0(outdir, args[1], "_harmonised_corrected.rds"))