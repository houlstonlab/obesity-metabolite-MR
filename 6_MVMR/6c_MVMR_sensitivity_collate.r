###
# Harmonise metabolite and obesity measure GWAS summary statistics and format
# for metaCCA analysis.
#
# args[1] = obesity measure
# args[2, 3, ...] = metabolite ids
###

library(tidyverse)
library(TwoSampleMR)

args <- commandArgs(trailingOnly = TRUE)

outdir <- "/path/to/your/output/directory/"
exposure <- args[1]
mediators <- args[-1]

# find which SNPs are present and pass harmonisation in all mediator GWAS
for (i in seq_along(mediators)) {
    # load list of all SNPs for a particular mediator
    exp_snps <- read_delim(paste0("/path/to/metabolite/exposures/", mediators[i], "/summary_stats.tsv"))

    # add exposure column
    exp_snps$exposure <- mediators[i]
    mediator_dat <- exp_snps %>% select(rsid, effect_allele, other_allele, eaf, beta, se, exposure)

    # remove SNPs with eaf > 0.99 or eaf < 0.01
    mediator_dat <- mediator_dat %>% filter(eaf >= 0.01, eaf <= 0.99)

    if (i == 1) {
        # cast in TwoSampleMR format
        mediator_dat <- format_data(
            mediator_dat,
            type = "exposure",
            phenotype_col = "exposure",
            snp_col = "rsid",
            beta_col = "beta",
            se_col = "se",
            eaf_col = "eaf",
            effect_allele_col = "effect_allele",
            other_allele_col = "other_allele"
        )

        harmonised_mediator <- mediator_dat
    } else {
        # cast in TwoSampleMR format
        mediator_dat <- format_data(
            mediator_dat,
            type = "outcome",
            phenotype_col = "exposure",
            snp_col = "rsid",
            beta_col = "beta",
            se_col = "se",
            eaf_col = "eaf",
            effect_allele_col = "effect_allele",
            other_allele_col = "other_allele"
        )

        harmonised_mediator <- harmonise_data(harmonised_mediator, mediator_dat)
        harmonised_mediator <- select(harmonised_mediator, c("SNP",
            "effect_allele.exposure", "other_allele.exposure", "beta.exposure",
            "se.exposure", "eaf.exposure", "exposure", "pval.exposure",
            "mr_keep", "pval_origin.exposure", "id.exposure"))
        harmonised_mediator <- filter(harmonised_mediator, mr_keep == TRUE)
    }
}

# find which SNPs are also present in the obesity measure GWAS
exposure_dat <- read_delim("/path/to/obesity/measure/MVMR_snps.tsv") %>%
                    mutate(exposure = "obesity_measure")
                    filter(eaf >= 0.01, eaf <= 0.99)

# cast in TwoSampleMR format
exposure_dat <- format_data(
    exposure_dat,
    type = "outcome",
    phenotype_col = "exposure",
    snp_col = "rsid",
    beta_col = "beta",
    se_col = "se",
    eaf_col = "eaf",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    pval_col = "pval"
)

harmonised_mediator <- harmonise_data(harmonised_mediator, exposure_dat)
harmonised_mediator <- select(harmonised_mediator, c("SNP",
    "effect_allele.exposure", "other_allele.exposure", "beta.exposure",
    "se.exposure", "eaf.exposure", "exposure", "pval.exposure",
    "mr_keep", "pval_origin.exposure", "id.exposure"))
harmonised_mediator <- filter(harmonised_mediator, mr_keep == TRUE)

# create a dataframe to collate harmonised data in metaCCA input format
# we add the first metabolite separately since it is the 'exposure' in the
# harmonised_mediator dataframe
collated_harmonised <- data.frame(allele_0 = rep(NA, nrow(harmonised_mediator)),
                                  allele_1 = rep(NA, nrow(harmonised_mediator)),
                                  V1 = rep(NA, nrow(harmonised_mediator)),
                                  V2 = rep(NA, nrow(harmonised_mediator)))
names(collated_harmonised)[names(collated_harmonised) == "V1"] <- paste0("trait", args[2], "_b")
names(collated_harmonised)[names(collated_harmonised) == "V2"] <- paste0("trait", args[2], "_se")

collated_harmonised$allele_0 <- harmonised_mediator$effect_allele.exposure
collated_harmonised$allele_1 <- harmonised_mediator$other_allele.exposure
collated_harmonised[[paste0("trait", args[2], "_b")]] <- harmonised_mediator$beta.exposure
collated_harmonised[[paste0("trait", args[2], "_se")]] <- harmonised_mediator$se.exposure

rownames(collated_harmonised) <- harmonised_mediator$SNP

# add the rest of the metabolites to the collated_harmonised dataframe
for (i in seq_along(mediators)) {
    if (i == 1) next

    # load list of all SNPs for a particular mediator
    exp_snps <- read_delim(paste0("/path/to/metabolite/exposures/", mediators[i], "/summary_stats.tsv"))

    # add exposure column
    exp_snps$exposure <- mediators[i]
    exp_snps <- exp_snps %>% select(rsid, effect_allele, other_allele, beta, se, exposure)

    # keep only SNPs that are in harmonised_mediator
    exp_snps <- filter(exp_snps, rsid %in% rownames(collated_harmonised))
    exp_snps <- filter(exp_snps, effect_allele == collated_harmonised[rsid, ]$allele_0 |
                            effect_allele == collated_harmonised[rsid, ]$allele_1)

    # this should remove the same duplicate SNPs as above
    mediator <- format_data(
        exp_snps,
        type = "exposure",
        phenotype_col = "exposure",
        snp_col = "rsid",
        beta_col = "beta",
        se_col = "se",
        effect_allele_col = "effect_allele",
        other_allele_col = "other_allele"
    )

    # convert back to tibble to allow use of mutate function
    mediator <- tibble(mediator)

    # if the effect and other allele are opposite, flip the beta
    mediator <- mutate(mediator, beta.exposure = ifelse(effect_allele.exposure == collated_harmonised[SNP, ]$allele_0 |
        effect_allele.exposure == collated_harmonised[SNP, ]$allele_1,
        ifelse(effect_allele.exposure == collated_harmonised[SNP, ]$allele_0, beta.exposure, -beta.exposure),
        stop("Effect allele does not match either allele in harmonised data.")))

    # format to match metaCCA input
    mediator <- mediator %>% arrange(match(SNP, rownames(collated_harmonised)))
    collated_harmonised[[paste0("trait", mediators[i], "_b")]] <- mediator$beta.exposure
    collated_harmonised[[paste0("trait", mediators[i], "_se")]] <- mediator$se.exposure
}

# add the obesity measure to the collated_harmonised dataframe
exposure_dat <- read_delim("/path/to/obesity/measure/MVMR_snps.tsv") %>%
                    mutate(exposure = "obesity_measure") %>%
                    select(rsid, effect_allele, other_allele, beta, se, exposure)

# keep only SNPs that are in harmonised_mediator
exposure_dat <- filter(exposure_dat, rsid %in% rownames(collated_harmonised))
exposure_dat <- filter(exposure_dat, effect_alelle == collated_harmonised[SNP, ]$allele_0 |
                        effect_allele == collated_harmonised[SNP, ]$allele_1)

# cast in TwoSampleMR format
exposure_dat <- format_data(
        exposure_dat,
        type = "exposure",
        phenotype_col = "exposure",
        snp_col = "rsid",
        beta_col = "beta",
        se_col = "se",
        effect_allele_col = "effect_allele",
        other_allele_col = "other_allele"
)

# convert back to tibble to allow use of mutate function
exposure_dat <- tibble(exposure_dat)

# if the effect and other allele are opposite, flip the beta
exposure_dat <- mutate(exposure_dat, beta.exposure = ifelse(effect_allele.exposure == collated_harmonised[SNP, ]$allele_0 |
    effect_allele.exposure == collated_harmonised[SNP, ]$allele_1,
    ifelse(effect_allele.exposure == collated_harmonised[SNP, ]$allele_0, beta.exposure, -beta.exposure),
    stop("Effect allele does not match either allele in harmonised data.")))

# format to match metaCCA input
exposure_dat <- exposure_dat %>% arrange(match(SNP, rownames(collated_harmonised)))
collated_harmonised[[paste0("trait", exposure, "_b")]] <- exposure_dat$beta.exposure
collated_harmonised[[paste0("trait", exposure, "_se")]] <- exposure_dat$se.exposure

write.table(collated_harmonised, paste0(outdir, "/collated_harmonised.tsv"), sep = "\t")