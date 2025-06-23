###
# Reclump union of exposure/mediator IVs and harmonise with outcome data for use in
# multivariable MR analysis.
#
# args[1, 2, 3, ...] = metabolite ids
###

library(tidyverse)
library(TwoSampleMR)

args <- commandArgs(trailingOnly = TRUE)

cancer <- "CRC"
outdir <- "/path/to/your/output/directory/"

# create vector of pruned GWAS files
filenames_exposure <- character(length(args) + 1)
for (i in seq_along(args)) {
    filenames_exposure[i] <- paste0("path/to/metabolite/exposures/", args[i], "/MVMR_snps.tsv")
}
filenames_exposure[length(args) + 1] <- "/path/to/obesity/measure/MVMR_snps.tsv"

# reclump union of all exposures' IVs
exposure <- mv_extract_exposures_local(
    filenames_exposure,
    sep = "\t",
    phenotype_col = "exposure",
    snp_col = "rsid",
    beta_col = "beta",
    se_col = "se",
    eaf_col = "eaf",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    pval_col = "pval",
    pval_threshold = 5e-08,
    plink_bin = "/path/to/plink/binary",
    bfile = "/path/to/10kg/reference",
    clump_r2 = 0.01,
    clump_kb = 500
)

# load cancer GWAS summary statistics
outcome <- read_delim("/path/to/CRC/GWAS/summary_stats.tsv") %>%
                arrange(pval) %>%
                filter(eaf >= 0.01, eaf <= 0.99)

outcome <- format_data(
    outcome,
    type = "outcome",
    snp_col = "rsid",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    pval_col = "pval",
    beta_col = "beta",
    se_col = "se",
    eaf_col = "eaf",
    samplesize_col = "N"
)

# harmonise exposure and outcome data and save to file
harmonised <- mv_harmonise_data(exposure, outcome)
saveRDS(harmonised, file = paste0(outdir, cancer, "_harmonised.rds"))