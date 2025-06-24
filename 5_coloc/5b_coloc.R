###
# Script to perform coloc analysis using metabolite and cancer GWAS.
#
# args[1] = path to metabolite GWAS folder
# args[2] = metabolite id
# args[3] = name of GWAS summary stats file
# args[4] = path to cancer file
# args[5] = number of cancer cases
# args[6] = number of cancer controls
###

library(tidyverse)
library(TwoSampleMR)
library(coloc)

args <- commandArgs(trailingOnly = TRUE)

exp_dir <- args[1]
exp <- args[2]
gwas_file <- args[3]
cancer_file <- args[4]
ncase <- as.numeric(args[5])
ncontrol <- as.numeric(args[6])

# load and format cancer outcome data
outcome <- read_delim(cancer_file)
outcome <- outcome %>% filter(eaf >= 0.01, eaf <= 0.99)
outcome <- format_data(
    outcome,
    type = "outcome",
    snp_col = "rsid",
    beta_col = "beta",
    se_col = "se",
    eaf_col = "eaf",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    pval_col = "pval",
    samplesize_col = "N",
    pos_col = "pos",
    chr_col = "chr"
)

# load SNPs acting as independent variables for the exposure
snp_list <- read_delim(paste0(exp_dir, "/coloc_snp.tsv"))

# load LD for SNPs within 500 kb of each independent SNP
ld <- read.table(paste0(exp_dir, "/r2_combined.ld"), header = TRUE)

for (chr in 1:22) {
    if (nrow(snp_list[snp_list$chr == chr, ]) == 0) next # skip if no IVs for this chromosome

    # load exposure GWAS summary stats for the chromosome
    exposure <- read_delim(paste0(exp_dir, "/", gwas_file, chr, ".tsv")) %>%
                    filter(eaf >= 0.01, eaf <= 0.99)

    exposure <- format_data(
        exposure,
        type = "exposure",
        snp_col = "rsid",
        beta_col = "beta",
        se_col = "se",
        eaf_col = "eaf",
        effect_allele_col = "effect_allele",
        other_allele_col = "other_allele",
        pval_col = "pval",
        samplesize_col = "N",
        pos_col = "pos",
        chr_col = "chr"
    )

    for (ind_snp in split(snp_list[snp_list$chr == chr, ], seq_len(nrow(snp_list[snp_list$chr == chr, ])))) {

        # filter for SNPs within 500 kb of independent SNP
        exposure_filtered <- exposure %>%
            filter(between(pos.exposure, ind_snp$pos - 500000, ind_snp$pos + 500000)
            & chr.exposure == ind_snp$chr)

        # harmonise exposure and outcome data
        harmonised <- harmonise_data(exposure_filtered, outcome)
        print(paste0("Number of SNPs in harmonised data ", nrow(harmonised)))

        # calculate prior probabilities according to rules given on
        # https://chr1swallace.shinyapps.io/coloc-priors/
        p1 <- 1 / (10 * nrow(harmonised))
        p2 <- p1
        p12 <- p1 / 10

        # create coloc datasets and check them
        coloc_exp <- list(beta = harmonised$beta.exposure, varbeta = harmonised$se.exposure^2,
            snp = harmonised$SNP, position = harmonised$pos.exposure, N = harmonised$samplesize.exposure,
            pvalues = harmonised$pval.exposure)
        coloc_exp$type <- "quant" # metabolites are continuous traits
        coloc_exp$sdY <- 1 # metabolite GWAS have been standardised to have a standard deviation of 1

        coloc_out <- list(beta = harmonised$beta.outcome, varbeta = harmonised$se.outcome^2,
            snp = harmonised$SNP, position = harmonised$pos.outcome, N = harmonised$samplesize.outcome,
            pvalues = harmonised$pval.outcome)
        coloc_out$type <- "cc" # cancer is a case-control trait

        check_dataset(coloc_exp)
        check_dataset(coloc_out)

        # perform coloc analysis
        res <- coloc.abf(dataset1 = coloc_exp, dataset2 = coloc_out,
                        p1 = p1, p2 = p2, p12 = p12)
        print(res)
        saveRDS(res, paste0(exp_dir, "/", exp, "_", ind_snp$rsid, "_res.rds"))

        print("Saving top 10 probable causal SNPs to file")
        res$results %>%
            slice_max(SNP.PP.H4, n = 10) %>%
            write_tsv(paste0(exp_dir, "/", exp, "_", ind_snp$rsid, "_coloc.tsv"))

        print("Results for H4 > 0.01")
        print(subset(res$results, SNP.PP.H4 > 0.01))

        print("95% credible set")
        o <- order(res$results$SNP.PP.H4, decreasing = TRUE)
        cs <- cumsum(res$results$SNP.PP.H4[o])
        w <- which(cs > 0.95)[1]
        print(res$results[o, ][1:w, ]$snp)

        print("Sensitivity analysis")
        pdf(paste0(exp_dir, "/sens_analysis_", ind_snp$rsid, ".pdf"))
        sensitivity(res, rule = "H4 > 0.8")
        dev.off()
    }
}
