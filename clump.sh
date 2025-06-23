#!/bin/bash

###
# Clump GWAS significant SNPs for a metabolite exposure to find independent variables.
###

exp="metabolite_id"

# create file only containing GWAS significant (p<5e-8) SNPs
Rscript find_gwas_sig_snps.R /path/to/metabolite/exposures/${exp}/summary_stats.tsv ${exp}

for chr in $(seq 1 22)
do

# clump GWAS significant SNPs to find independent SNPs
plink --bfile /path/to/10kg/reference/1kg_uk10k_eur_dedup/1kg_uk10k_eur_chr${chr}_filtered_dedup \
      --clump chr${chr}_gwas_sig_snps.txt \
      --clump-p1 5e-8 \
      --clump-p2 5e-8 \
      --clump-r2 0.01 \
      --clump-kb 500 \
      --clump-field pval \
      --clump-snp-field rsid \
      --clump-allow-overlap \
      --noweb \
      --out ${exp}_chr${chr}

done