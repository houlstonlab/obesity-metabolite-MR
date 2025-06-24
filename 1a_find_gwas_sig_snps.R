###
# This script filters GWAS summary statistics for a given metabolite to find significant SNPs.
#
# args[1] = path to metabolite GWAS summary stats
# args[2] = metabolite id
###

library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

# read in GWAS summary stats
data <- read_delim(args[1])
exp <- args[2]

for (chr in 1:22) {
    # filter SNPs to keep only those with:
    # 1. pval < 5e-8
    # 2. maf between 0.01 and 0.99
    # 3. not in HLA region
    data <- data %>%
        filter(chr == chr) %>%
        filter(pval < 5e-8) %>%
        filter(eaf >= 0.01, eaf <= 0.99) %>%
        filter(!(between(pos, 28477897, 33448354) & chr == 6)) # assumes GRCh37 format
    write_tsv(data, paste0("/path/to/metabolite/exposures/", args[2], "/chr", chr, "_gwas_sig_snps.txt"))
}