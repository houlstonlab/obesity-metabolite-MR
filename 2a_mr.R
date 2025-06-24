###
# Perform two-sample Mendelian Randomization (MR) analysis using
# metabolite summary GWAS data.
#
# The script may need to be edited to suit different exposure/outcome datasets.
#
# args[1] = cancer name
# args[2] = lifetime risk of particular cancer
# args[3] = ncase
# args[4] = ncontrol
# args[5] = outcome data
# args[6] = outcome pval column name
###

library(tidyverse)
library(TwoSampleMR)

args <- commandArgs(trailingOnly = TRUE)

### section to edit ###
# assign command line arguments to variables
ncase <- as.integer(args[3])
ncontrol <- as.integer(args[4])
lifetimerisk <- as.numeric(args[2])
cancer <- args[1]
outdir <- "/path/to/your/output/directory/"

# find all metabolite exposures
listexp <- list.dirs(
    "/path/to/metabolite/exposures/",
    recursive = FALSE,
    full.names = FALSE
)

# create empty tibble to store exposure IVs
exposure <- tibble(
    chr = integer(),
    pos = integer(),
    rsid = character(),
    effect_allele = character(),
    other_allele = character(),
    N = integer(),
    eaf = numeric(),
    beta = numeric(),
    se = numeric(),
    pval = numeric(),
    exposure = character()
)

# collate IVs for each exposure into a single tibble
for (dir in listexp) {
    for (chr in 1:22) {
        # PLINK does not create a .clumped file for exposures without any significant SNPs
        if (!file.exists(paste0("/path/to/metabolite/exposures/", dir, "/", dir, "_chr", chr, ".clumped"))) next

        # load list of all GWAS significant SNPs for a particular exposure and chromosome
        exp_chr_sig_snps <- read_delim(paste0("/path/to/metabolite/exposures/", dir, "/chr", chr, "_gwas_sig_snps.txt"),
                                col_types = "iiicccidddddc")
        # load list of lead SNPs for a particular exposure and chromosome
        exp_chr_clump <- read_table(paste0("/path/to/metabolite/exposures/", dir, "/", dir, "_chr", chr, ".clumped"),
                                col_types = "iicidiiiiiic")

        # drop any IVs without rsid or position
        exp_chr_clump <- exp_chr_clump %>% drop_na()
        exp_chr_sig_snps <- exp_chr_sig_snps %>% filter(paste(rsid, pos) %in% paste(exp_chr_clump$SNP, exp_chr_clump$BP))

        # add exposure IVs to exposure tibble
        if (nrow(exp_chr_sig_snps) > 0) {
            # attach the name of the exposure to the IVs
            exp_chr_sig_snps$exposure <- dir
            exposure <- bind_rows(exposure, exp_chr_sig_snps)
        }
    }
}

# order rows by ascending pval so that lowest pval is kept by TwoSampleMR in
# the case of duplicated SNPs
exposure <- exposure %>% arrange(pval)

# remove SNPs with eaf > 0.99 or eaf < 0.01
# shouldn't be needed as this is performed in the clumping step, but good to check
exposure <- exposure %>% filter(eaf >= 0.01, eaf <= 0.99)

# cast in TwoSampleMR format
exposure <- format_data(
        exposure,
        type = "exposure",
        phenotype_col = "exposure",
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

# read in outcome data
outcome <- read_delim(args[5])

# order rows by ascending pval so that lowest pval is kept by TwoSampleMR in
# the case of duplicated SNPs
outcome <- outcome %>% arrange(args[6])

# remove SNPs with eaf > 0.99 or eaf < 0.01
outcome <- outcome %>% filter(eaf >= 0.01, eaf <= 0.99)

outcome <- format_data(
    outcome,
    type = "outcome",
    phenotype_col = "exposure",
    snp_col = "rsid",
    beta_col = "beta",
    se_col = "se",
    eaf_col = "eaf",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    pval_col = args[6],
    samplesize_col = "N",
    pos_col = "pos",
    chr_col = "chr"
)

### end of section to edit ###

# harmonise the exposure and outcome data, check it looks sensible, then write table
harmonised <- harmonise_data(exposure, outcome)

head(harmonised)
write_tsv(harmonised, paste0(outdir, cancer, "_harmonised.tsv"))


# we need a minor allele frequency (maf) column for f-statistic and power calculations.
# sometimes the effect column isn't always the minor allele
harmonised <- harmonised %>%
    mutate(
        maf = if_else(eaf.exposure > 0.5, 1 - eaf.exposure, eaf.exposure)
	)


# compute each SNPs PVE and add on column
# use PVE equation from http://journals.plos.org/plosone/article/file?type=supplementary&id=info:doi/10.1371/journal.pone.0120758.s001
harmonised <-  harmonised %>%
    mutate(pve =
		(2 * beta.exposure ^ 2 * maf * (1 - maf)) /
        ((2 * beta.exposure ^ 2 * maf * (1 - maf)) +
        (se.exposure ^ 2 * 2 * samplesize.exposure * maf * (1 - maf)))
	)


# compute f-statistic for a trait
# uses F statistic equations from https://academic.oup.com/ije/article/40/3/740/741448/Power-and-instrument-strength-requirements-for
exposure_statistics <- harmonised %>%
	group_by(id.exposure) %>%
	summarise(
		exposure = first(exposure),
		n_snp = n(),
		r2 = sum(pve),
		samplesize = first(samplesize.exposure)
	) %>%
	mutate(fstat = (r2 * (samplesize - 1 - n_snp)) / ((1 - r2) * n_snp))


# calculate power for a range of odds ratios
# uses power equations from https://cnsgenomics.shinyapps.io/mRnd/
# first make function for power
ors <- c(1.05, 1.1, 1.25, 1.33, 1.5)
n <- ncase + ncontrol
k <- ncase / n
alpha <- 0.05
z_alpha <- qnorm(1 - alpha / 2, 0, 1)


calculate_power <-  function(r2, n, k, ors) {
	b01 <- k * (ors / (1 + k * (ors - 1)) - 1)
    ncp <- (n * r2) * (k * (ors / (1 + k * (ors - 1)) - 1)) ^ 2 /
		(k * (1 - k) - b01 ^ 2)
    power <- 1 + pnorm(-z_alpha - sqrt(ncp), 0, 1) -
		pnorm(z_alpha - sqrt(ncp), 0, 1)
    power
	}

# calculate power for each odds ratio and bind to exposure statistics
power_table <- map_dfc(
	ors, ~mutate(exposure_statistics, !!sym(paste0("power", .)) :=
	calculate_power(r2, n, k, .)) %>%
	select(starts_with("power"))
)

# bind fstat and pve etc to power table and check it looks sensible
stats_dat <- bind_cols(exposure_statistics, power_table)
head(stats_dat)
dim(stats_dat)

# write out table
write_tsv(stats_dat, paste0(outdir, "power_table_", cancer, ".txt"))

# filter to remove NAs
stats_dat_no_na <- stats_dat %>% filter(!is.na(power1.05))

# filter to remove SNPs with fstat < 10
filter_data_f10 <- stats_dat_no_na %>% filter(fstat >= 10)

# filter to remove duplicate exposures, taking one with highest power
filter_data_f10 %>%
	group_by(exposure) %>%
	arrange(desc(power1.05)) %>%
	filter_data_maxp <- slice(1)
dim(filter_data_maxp)

write_tsv(filter_data_maxp, paste0(outdir, "power_table_", cancer, "_filter.txt"))

power_id <- filter_data_maxp %>%
            ungroup %>%
            select(1)

pruned_harmonised_data <- left_join(power_id, harmonised, by = "id.exposure")
write_tsv(pruned_harmonised_data, paste0(outdir, cancer, "_filtered_harmonised_data.txt"))


# run mr
mr_methods <- c(
	"mr_wald_ratio", "mr_ivw", "mr_ivw_fe", "mr_two_sample_ml",
	"mr_simple_median", "mr_weighted_median", "mr_simple_mode",
	"mr_weighted_mode", "mr_egger_regression"
)
mr_res <- mr(pruned_harmonised_data, method_list = mr_methods)
write_tsv(mr_res, paste0(outdir, cancer, "_mr_results.txt"))


# run sensitivity tests
het <- mr_heterogeneity(pruned_harmonised_data)
write_tsv(het, paste0(outdir, cancer, "_heterogeneity_MR.txt"))

pleio <- mr_pleiotropy_test(pruned_harmonised_data)
write_tsv(pleio, paste0(outdir, cancer, "_pleiotropy_egger_MR.txt"))

loo <- mr_leaveoneout(pruned_harmonised_data)
write_tsv(loo, paste0(outdir, cancer, "_leaveoneout_MR.txt"))

single <- mr_singlesnp(pruned_harmonised_data)
write_tsv(single, paste0(outdir, cancer, "_single_snp_MR.txt"))

#these are the values which decide what will be flagged in the analysis
heterogeneity_flag_pval <- 0.05
snp_flag_min <- 2
egger_flag_pval <- 0.05
leave_one_out_flag_pval <- 0.05

# flag exposures based on sensitivity analysis results
exposures_flags <- list()
exposures_flags$few_snps <- mr_res$exposure[mr_res$nsnp < snp_flag_min]
# ifelse statement hopefully ensures column headers are written even when
# the particular test could not be run
exposures_flags$pleio <-
	ifelse(pleio$exposure[pleio$pval < egger_flag_pval & !is.na(pleio$pval)] != NULL,
	pleio$exposure[pleio$pval < egger_flag_pval & !is.na(pleio$pval)], NA)
exposures_flags$heterogeneity <-
    ifelse(het$exposure[het$method == "Inverse variance weighted" & het$Q_pval <
	heterogeneity_flag_pval & !is.na(het$Q_pval)] != NULL,
	het$exposure[het$method == "Inverse variance weighted" & het$Q_pval <
	heterogeneity_flag_pval & !is.na(het$Q_pval)], NA)

# function to identify exposures to flag due to leave-one-out sensitivity analysis
# if pval < p when all SNPs considered, but the removal of any single results in pval > p, then flag
identify.leave.one.out.flag <- function(res, p) {
	res_split <- split(res, res$exposure)
    names(res_split)[sapply(res_split,
		function(res_exposure) {
			res_exposure$p[res_exposure$SNP == "All"] < p &
			!all(res_exposure$p[res_exposure$SNP != "All"] < p)
		})]
}

exposures_flags$leave_one_out_ivw <- identify.leave.one.out.flag(loo, leave_one_out_flag_pval)

# join flags to mr results and write out table
for (flag in names(exposures_flags)) mr_res[[paste("flag", flag, sep = ".")]] <-
	mr_res$exposure %in% exposures_flags[[flag]]
# change value of flags for exposure with nsnp < 2 to NA from FALSE as these
# tests cannot be run with fewer than 2 SNPs
mr_res[["flag.pleio"]][mr_res$nsnp < 2] <- NA
mr_res[["flag.heterogeneity"]][mr_res$nsnp < 2] <- NA
write_tsv(mr_res, paste0(outdir, cancer, "_flagged_results_mr_all.txt"))