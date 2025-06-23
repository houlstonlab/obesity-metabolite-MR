#!/bin/bash

###
# Example batch script to run Mendelian Randomization (MR) analyses for multiple cancers using metabolomics data.
###

home_dir="/path/to/home/directory"

cd ${home_dir}

# list of cancers
declare -a cancers=("CRC" "RCC" "breast" "endometrial" "lung" "oesophageal" "ovarian" "prostate")
cancerslength=${#cancers[@]}

# list of lifetime risk
declare -a liferisks=("0.06" "0.02" "0.13" "0.03" "0.06" "0.015" "0.02" "0.13")

# list of ncases
declare -a ncases=("73673" "10784" "133384" "8758" "29266" "16790" "26293" "79194")

# list of ncontrols
declare -a ncontrols=("86854" "20406" "113789" "46126" "56450" "32476" "68502" "61112")

# list of cancer summary stats
declare -a summarys=(
    "/path/to/CRC/GWAS/summary_stats.tsv"
    "/path/to/RCC/GWAS/summary_stats.tsv"
    "/path/to/breast/GWAS/summary_stats.tsv"
    "/path/to/endometrial/GWAS/summary_stats.tsv"
    "/path/to/lung/GWAS/summary_stats.tsv"
    "/path/to/oesophageal/GWAS/summary_stats.tsv"
    "/path/to/ovarian/GWAS/summary_stats.tsv"
    "/path/to/prostate/GWAS/summary_stats.tsv")

# list of pval column headers
declare -a pvals=("P_value" "P_value" "p.meta" "pval" "p_value" "P-value" "pval.outcome" "Pvalue")

# loop over cancers
for (( i=0; i<${cancerslength}; i++ ))
do

if [ ! -d ${cancers[$i]} ]; then
  mkdir ${cancers[$i]}
fi

cd ${cancers[$i]}

cp ${home_dir}/mr.R ./

cat << EOF > mr.sh
#!/bin/bash

#SBATCH --job-name=MR_metabolomics_${cancers[$i]}
#SBATCH --partition=compute
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=4:00:00
#SBATCH --output=${cancers[$i]}%j.o
#SBATCH --error=${cancers[$i]}%j.e

Rscript ${home_dir}/${cancers[$i]}/mr.R ${cancers[$i]} ${liferisks[$i]} ${ncases[$i]} ${ncontrols[$i]} ${summarys[$i]} ${pvals[$i]}
EOF

sbatch mr.sh

cd ${home_dir}

done