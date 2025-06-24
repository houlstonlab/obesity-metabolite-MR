###
# Collate independent variables for a metabolite exposure for coloc analysis.
#
# args[1] = path to metabolite GWAS folder
# args[2] = metabolite id
###

library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

metabolite_folder <- args[1]
exp <- args[2]

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

for (chr in 1:22) {
    # PLINK does not create a .clumped file for exposures without any significant SNPs
    if (!file.exists(paste0(metabolite_folder, exp, "/", exp, "_chr", chr, ".clumped"))) next

    # load list of all GWAS significant SNPs for a particular exposure and chromosome
    exp_chr_sig_snps <- read_delim(paste0(metabolite_folder, exp, "/chr", chr, "_gwas_sig_snps.txt"),
                            col_types = "iiicccidddddc")
    # load list of lead SNPs for a particular exposure and chromosome
    exp_chr_clump <- read_table(paste0(metabolite_folder, exp, "/", exp, "_chr", chr, ".clumped"),
                            col_types = "iicidiiiiiic")

    # drop any IVs without rsid or position
    exp_chr_clump <- exp_chr_clump %>% drop_na()
    exp_chr_sig_snps <- exp_chr_sig_snps %>% filter(paste(rsid, pos) %in% paste(exp_chr_clump$SNP, exp_chr_clump$BP))

    # add exposure IVs to clump_dat tibble
    if (nrow(exp_chr_sig_snps) > 0) {
        clump_dat <- bind_rows(clump_dat, exp_chr_sig_snps)
    }
}

write(clump_dat$SNP, paste0(metabolite_folder, exp, "/", "plink_ind_snp.txt"))
write_tsv(clump_dat, paste0(metabolite_folder, exp, "/", "coloc_snp.tsv"))
