###
# Post-hoc power calculations for quantitative outcomes.
###

library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

cancer <- args[1]
outdir <- "/path/to/your/output/directory/"

harmonised <- read_delim(paste0(outdir, cancer, "_harmonised.tsv"))
results <- read_delim(paste0(outdir, cancer, "_flagged_results_mr_all.txt"))

# we need a minor allele frequency (maf) column for power calculations.
# sometimes the effect column isn't always the minor allele
harmonised <- harmonised %>%
    mutate(
        maf = if_else(eaf.exposure > 0.5, 1 - eaf.exposure, eaf.exposure)
	)

## compute each SNPs PVE and total r2 and variance
## use PVE equation from http://journals.plos.org/plosone/article/file?type=supplementary&id=info:doi/10.1371/journal.pone.0120758.s001
harmonised <- harmonised %>%
    mutate(pve.exposure =
		(2 * beta.exposure ^ 2 * maf * (1 - maf)) /
        ((2 * beta.exposure ^ 2 * maf * (1 - maf)) +
        (se.exposure ^ 2 * 2 * samplesize.exposure * maf * (1 - maf)))
	) %>%
    mutate(pve.outcome = 
        (2 * beta.outcome ^ 2 * maf * (1 - maf)) /
        ((2 * beta.outcome ^ 2 * maf * (1 - maf)) +
        (se.outcome ^ 2 * 2 * samplesize.outcome * maf * (1 - maf)))
    )

alpha <- 0.05
z_alpha <- qnorm(1 - alpha / 2, 0, 1)

# power for continous outcome calculated using https://academic.oup.com/ije/article/42/5/1497/623616
exposure_statistics <- harmonised %>%
	group_by(id.outcome, id.exposure) %>%
	summarise(
		exposure = first(exposure),
        outcome = first(outcome),
		n_snp = n(),
		r2.exposure = sum(pve.exposure),
        r2.outcome = sum(pve.outcome),
		samplesize = first(samplesize.exposure),
        var_x = sum((2 * beta.exposure ^ 2 * maf * (1 - maf)) / r2.exposure),
        var_y = sum((2 * beta.outcome ^ 2 * maf * (1 - maf)) / r2.outcome),
        beta = first(results$b[(results$exposure == exposure) & (results$outcome == outcome)]), # should always be IVW-RE or Wald ratio
        #beta = 0.049, OR = 1.05; beta = 0.095, OR = 1.1; beta = 0.22, OR = 1.25; beta = 0.26, OR = 1.33; beta = 0.41, OR = 1.5
        ncp_0.049 = (samplesize * r2.exposure * var_x) * (0.049 + (beta - 0.049) * var_x / # 0.049/0.095/0.22/etc. here are the effect sizes, not odds ratios
            (samplesize * r2.exposure)) ^ 2 / (var_y - var_x * 0.049 * (2 * beta - 0.049)),
        ncp_0.095 = (samplesize * r2.exposure * var_x) * (0.095 + (beta - 0.095) * var_x / 
            (samplesize * r2.exposure)) ^ 2 / (var_y - var_x * 0.095 * (2 * beta - 0.095)),
        ncp_0.22 = (samplesize * r2.exposure * var_x) * (0.22 + (beta - 0.22) * var_x / 
            (samplesize * r2.exposure)) ^ 2 / (var_y - var_x * 0.22 * (2 * beta - 0.22)),
        ncp_0.26 = (samplesize * r2.exposure * var_x) * (0.26 + (beta - 0.26) * var_x / 
            (samplesize * r2.exposure)) ^ 2 / (var_y - var_x * 0.26 * (2 * beta - 0.26)),
        ncp_0.41 = (samplesize * r2.exposure * var_x) * (0.41 + (beta - 0.41) * var_x / 
            (samplesize * r2.exposure)) ^ 2 / (var_y - var_x * 0.41 * (2 * beta - 0.41)),
        power_0.049 = 1 + pnorm(-z_alpha - sqrt(ncp_0.049), 0, 1) -
            pnorm(z_alpha - sqrt(ncp_0.049), 0, 1),
        power_0.095 = 1 + pnorm(-z_alpha - sqrt(ncp_0.095), 0, 1) -
            pnorm(z_alpha - sqrt(ncp_0.095), 0, 1),
        power_0.22 = 1 + pnorm(-z_alpha - sqrt(ncp_0.22), 0, 1) -
            pnorm(z_alpha - sqrt(ncp_0.22), 0, 1),
        power_0.26 = 1 + pnorm(-z_alpha - sqrt(ncp_0.26), 0, 1) -
            pnorm(z_alpha - sqrt(ncp_0.26), 0, 1),
        power_0.41 = 1 + pnorm(-z_alpha - sqrt(ncp_0.41), 0, 1) -
            pnorm(z_alpha - sqrt(ncp_0.41), 0, 1),
        .groups = "drop"
	)


write_tsv(exposure_statistics, paste0(outdir, "posthoc_quantitative_outcome_power.tsv"))
