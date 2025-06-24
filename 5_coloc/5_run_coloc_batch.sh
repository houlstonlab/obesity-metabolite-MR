#!/bin/bash

###
# Example batch script to run coloc analysis for metabolite exposures and cancer GWAS summary stats.
###

exp="metabolite_id"
cancer="CRC"
ncases=73673
ncontrols=86854

# calculate LD r^2 values around 500 kb of each independent SNP
cd /path/to/metabolite/exposures/${exp}/

# create list of independent SNPs and corresponding GWAS data
Rscript 5a_find_ind_snp.R /path/to/metabolite/exposures/ ${exp}

# calculate LD r^2 values for each independent SNP
for chr in $(seq 1 22)
do

plink --bfile /path/to/10kg/reference/1kg_uk10k_eur_dedup/1kg_uk10k_eur_chr${chr}_filtered_dedup \
      --r2 
      --ld-window 10000000 # set arbitrarily high to ensure all pairs are included
      --ld-window-kb 500 
      --ld-window-r2 0.01 
      --ld-snp-list /path/to/metabolite/exposures/${exp}/plink_ind_snp.txt 
      --out r2_chr${chr}

# combine all per chromosome LD files
rm -rf r2_combined.ld
echo $'CHR_A\tBP_A\tSNP_A\tCHR_B\tBP_B\tSNP_B\tR2' > r2_combined.ld
cat r2_chr${chr}.ld | awk 'NR>1' >> r2_combined.ld

done

Rscript 5b_coloc.R /path/to/metabolite/exposures/${exp} \
                ${exp} \
                summary_stats_chr \
                /path/to/${cancer}/GWAS/summary_stats.tsv ${ncases} ${ncontrols}