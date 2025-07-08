###
# Find proxy SNPs for missing exposure IVs in cancer GWAS data
#
# The script may need to be edited to suit different exposure/outcome datasets.
#
# args[1] = cancer name
# args[2] = cancer GWAS data file
###

library(tidyverse)
library(TwoSampleMR)

args <- commandArgs(trailingOnly = TRUE)

cancer <- args[1]
cancer_gwas <- args[2]
outdir <- "/path/to/your/output/directory/"

# combined exposure data for all metabolites
exposure <- read_delim("/path/to/metabolite/exposures/all_metabolite_ivs.txt")

# load harmonised data from MR without proxies
harmonised_data <- read_delim(paste0(outdir, cancer, "_harmonised.tsv"))

# find IVs that are not in the harmonised data
missing <- exposure[!exposure$SNP %in% harmonised_data$SNP, ]
write_tsv(missing, paste0(outdir, cancer, "_missing_exposure_snps.txt"))

# load cancer data
outcome <- read_delim(cancer_gwas)

# remove SNPs with eaf > 0.99 or eaf < 0.01
outcome <- outcome %>% filter(eaf >= 0.01, eaf <= 0.99)

# plink output format
proxy_snps <- tibble(
    CHR_A = integer(),
    BP_A = integer(),
    SNP_A = character(),
    MAF_A = numeric(),
    CHR_B = integer(),
    BP_B = integer(),
    SNP_B = character(),
    PHASE = character(),
    MAF_B = numeric(),
    R2 = numeric()
)

# cancer GWAS data format
proxy_snps_gwas <- tibble(
    rsid = character(),
    effect_allele = character(),
    other_allele = character(),
    pval = numeric(),
    beta = numeric(),
    se = numeric(),
    eaf = numeric()
)

for (chr in 1:22) {
    # find proxy SNPs per chromosome
    # if your reference data isn't split by chromosome, you can remove the loop
    write.table(unique(missing$SNP[missing$chr.exposure == chr]),
        paste0(outdir, cancer, "_missing_exposure_snps_chr", chr, ".txt"),
        row.names = FALSE, col.names = FALSE, quote = FALSE)
    system(paste0("plink --bfile /path/to/10kg/reference/1kg_uk10k_eur_dedup/1kg_uk10k_eur_chr",
        chr, "_filtered_dedup --r2 in-phase with-freqs gz --ld-window 100 ",
        "--ld-window-kb 500 --ld-window-r2 0.8 --ld-snp-list ", outdir, cancer,
        "_missing_exposure_snps_chr", chr, ".txt", " --out r2_chr", chr))
    
    # if no proxy SNPs found, skip to next chromosome
    if (!file.exists(paste0("r2_chr", chr, ".ld.gz"))) {
        next
    }

    # order proxy SNPs by r^2 and remove SNP_A == SNP_B
    r2_data <- read.table(paste0("r2_chr", chr, ".ld.gz"), header = TRUE)
    r2_data <- r2_data %>%
        filter(SNP_A != SNP_B) %>%
        arrange(SNP_A, desc(R2))

    # find proxy SNPs that are in the outcome data
    r2_data <- r2_data %>%
        filter(SNP_B %in% outcome$SNP)

    # keep the proxy SNP with highest r^2
    r2_data <- r2_data %>%
        group_by(SNP_A) %>%
        slice(1) %>%
        ungroup()
    proxy_snps <- bind_rows(proxy_snps, r2_data)

    # expand r2_data to allow for multiple SNP_A values per SNP_B
    proxy_map <- r2_data %>%
        select(SNP_A, SNP_B, PHASE)

    # PHASE is in the format SNP_A allele + SNP_B allele, however there is
    # no demarcation between the SNP_A and SNP_B allele
    # use .bim file to find SNP_A allele and work out SNP_B allele
    bim_file <- read_delim(paste0("/path/to/10kg/reference/1kg_uk10k_eur_dedup/1kg_uk10k_eur_chr",
        chr, "_filtered_dedup.bim"),
        col_names = c("chr", "SNP", "buffer", "pos", "a1", "a2")) %>%
	    separate_longer_delim(SNP, delim = ";")

    proxy_map <- left_join(proxy_map, bim_file, by = c("SNP_A" = "SNP"), keep = TRUE) %>%
	    select(-c("chr", "buffer", "pos", "SNP")) %>%
	    rename(SNP_A_a1 = a1, SNP_A_a2 = a2)

    proxy_map <- left_join(proxy_map, bim_file, by = c("SNP_B" = "SNP"), keep = TRUE) %>%
	    select(-c("chr", "buffer", "pos", "SNP")) %>%
	    rename(SNP_B_a1 = a1, SNP_B_a2 = a2)

    # if SNP_A_a1 + SNP_B_a1 != PHASE then need to swap SNP_B_a1 and SNP_B_a2
    proxy_map <- proxy_map %>% mutate(SNP_B_a1_copy = SNP_B_a1)
    proxy_map <- proxy_map %>% mutate(SNP_B_a1 = ifelse(paste0(SNP_A_a1, SNP_B_a1) %in% unlist(str_split(PHASE, "/")), SNP_B_a1, SNP_B_a2)) %>%
		 mutate(SNP_B_a2 = ifelse(paste0(SNP_A_a1, SNP_B_a1) %in% unlist(str_split(PHASE, "/")), SNP_B_a2, SNP_B_a1_copy)) %>%
		 select(-SNP_B_a1_copy)

    # join outcome data with proxy_map to allow multiple SNP_A per SNP_B
    proxy_snps_outcome <- outcome %>%
        filter(SNP %in% proxy_map$SNP_B) %>%
        select(rsid, effect_allele, other_allele, pval, beta, se, eaf) %>%
        left_join(proxy_map, by = c("rsid" = "SNP_B"), keep = TRUE) %>%
        mutate(proxy_SNP = rsid)

    # harmonise SNP_A and SNP_B
    proxy_snps_outcome <- proxy_snps_outcome %>% mutate(other_allele = ifelse(SNP_B_a1 == effect_allele, SNP_A_a2, SNP_A_a1)) %>%
	                      mutate(effect_allele = ifelse(SNP_B_a1 == effect_allele, SNP_A_a1, SNP_A_a2)) %>%
                          mutate(rsid = SNP_A) %>%
                          select(-c(SNP_A, SNP_B))

    # add to proxy_snps_gwas tibble
    if (nrow(proxy_snps_outcome) > 0) {
        proxy_snps_gwas <- bind_rows(proxy_snps_gwas, proxy_snps_outcome)
    }
}

# write proxy SNPs to file
write_tsv(proxy_snps, paste0(outdir, cancer, "_proxy_snps.tsv"))

# write proxy SNPs GWAS data to file
write_tsv(proxy_snps_gwas, paste0(outdir, cancer, "_proxy_snps_gwas.tsv"))
