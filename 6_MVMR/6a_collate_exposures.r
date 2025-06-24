###
# Used to shrink the metabolite GWAS data to only those SNPs that are IVs of any
# metabolite mediator or obesity measure.
# This vastly reduces the size of the data and speeds up the reclumping step.
#
# args[1, 2, 3, ...] = metabolite ids
###

library(tidyverse)
library(TwoSampleMR)

args <- commandArgs(trailingOnly = TRUE)

exposure_names <- args

# create tibble of all metabolite and obesity measure IVs
clump_dat <- tibble(
    chr = integer(),
    pos = integer(),
    rsid = character(),
    effect_allele = character(),
    other_allele = character(),
    N = integer(),
    eaf = numeric(),
    beta = numeric(),
    se = numeric(),
    pval = numeric()
)

# add IVs for each metabolite
for (exp_name in exposure_names) {
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
            clump_dat <- bind_rows(clump_dat, exp_chr_sig_snps)
        }
    }
}

# add IVs for obesity measure
for (chr in 1:22) {
    # PLINK does not create a .clumped file for exposures without any significant SNPs
    if (!file.exists(paste0("/path/to/obesity/measure/", dir, "/", dir, "_chr", chr, ".clumped"))) next

    # load list of all GWAS significant SNPs for a particular exposure and chromosome
    exp_chr_sig_snps <- read_delim(paste0("/path/to/obesity/measure/", dir, "/chr", chr, "_gwas_sig_snps.txt"),
                            col_types = "iiicccidddddc")
    # load list of lead SNPs for a particular exposure and chromosome
    exp_chr_clump <- read_table(paste0("/path/to/obesity/measure/", dir, "/", dir, "_chr", chr, ".clumped"),
                            col_types = "iicidiiiiiic")

    # drop any IVs without rsid or position
    exp_chr_clump <- exp_chr_clump %>% drop_na()
    exp_chr_sig_snps <- exp_chr_sig_snps %>% filter(paste(rsid, pos) %in% paste(exp_chr_clump$SNP, exp_chr_clump$BP))

    # add exposure IVs to exposure tibble
    if (nrow(exp_chr_sig_snps) > 0) {
        clump_dat <- bind_rows(clump_dat, exp_chr_sig_snps)
    }
}

# keep each IV only once
clump_dat <- clump_dat %>% distinct(rsid, .keep_all = TRUE)
clump_dat <- clump_dat %>% select(-c(N, eaf, beta, se, pval))

# next, we load each exposure one at a time and keep only SNPs in clump_dat
# this will include some SNPs in LD, however the clumping step of
# mv_extract_exposures_local will remove these

# create dataframe with flipped effect and non-effect alleles to account
# for SNPs that are not harmonised between clump_dat and raw exposure data
# could instead just search by rsid if you are sure that all SNPs in the raw
# exposure data have rsids
clump_dat_flipped <- rename(clump_dat, c(temp = "effect_allele",
                        effect_allele = "other_allele")) %>%
                        rename(other_allele = "temp")

# add SNPs in clump_dat to exposure tibble for each metabolite
for (exp_name in exposure_names) {
    exp_dat <- read_delim(paste0("/path/to/metabolite/exposures/", exp_name,
                                 "/summary_stats.tsv"), col_types = "iiicccidddddc")

    # filter by SNPs in clump_dat or clump_dat_flipped
    exp_dat <- exp_dat %>% filter(
        (paste(chr, pos, effect_allele, other_allele) %in%
         paste(clump_dat$chr, clump_dat$pos, clump_dat$effect_allele, clump_dat$other_allele)) |
        (paste(chr, pos, effect_allele, other_allele) %in%
         paste(clump_dat$chr, clump_dat$pos, clump_dat$other_allele, clump_dat$effect_allele)))

    # combine SNPs found in clump_dat and clump_dat_flipped
    exp_dat <- bind_rows(inner_join(exp_dat, clump_dat,
                         by = join_by(chr, pos, effect_allele, other_allele), keep = FALSE),
                         inner_join(exp_dat, clump_dat_flipped,
                         by = join_by(chr, pos, effect_allele, other_allele), keep = FALSE))

    # remove SNPs with eaf > 0.99 or eaf < 0.01
    exp_dat <- exp_dat %>% filter(eaf >= 0.01, eaf <= 0.99)

    exp_dat$exposure <- exp_name
    write_tsv(exp_dat, paste0("path/to/metabolite/exposures/", exp_name,
                              "/MVMR_snps.tsv"))
}

# add SNPs in clump_dat to exposure tibble for obesity measure
exp_dat <- read_delim("/path/to/obesity/measure/summary_stats.tsv", col_types = "iiicccidddddc")

# filter by SNPs in clump_dat
exp_dat <- exp_dat %>% filter(
    (paste(chr, pos, effect_allele, other_allele) %in%
     paste(clump_dat$chr, clump_dat$pos, clump_dat$effect_allele, clump_dat$other_allele)) |
    (paste(chr, pos, effect_allele, other_allele) %in%
     paste(clump_dat$chr, clump_dat$pos, clump_dat$other_allele, clump_dat$effect_allele)))

# combine SNPs found in clump_dat and clump_dat_flipped
exp_dat <- bind_rows(inner_join(exp_dat, clump_dat,
                     by = join_by(chr, pos, effect_allele, other_allele), keep = FALSE),
                     inner_join(exp_dat, clump_dat_flipped,
                     by = join_by(chr, pos, effect_allele, other_allele), keep = FALSE))

# remove SNPs with eaf > 0.99 or eaf < 0.01
exp_dat <- exp_dat %>% filter(eaf >= 0.01, eaf <= 0.99)

exp_dat$exposure <- "obesity_measure"
write_tsv(exp_dat, "/path/to/obesity/measure/MVMR_snps.tsv")
