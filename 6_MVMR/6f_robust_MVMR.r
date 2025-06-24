###
# Calculate direct effect size estimates for multiple exposures/mediators using
# IVW/Q-statistic minimsation MVMR.
###

library(metaCCA)
library(tidyverse)
library(MVMR)

outdir <- "/path/to/your/output/directory/"

S_XY <- read.table(paste0(outdir, "collated_harmonised.tsv"), header = TRUE, row.names = 1, stringsAsFactors = TRUE)
S_YY <- estimateSyy(S_XY = S_XY)

harmonised <- readRDS(paste0(outdir, "cancer_harmonised_corrected.rds"))

# save order of exposure/mediators
x_list <- colnames(harmonised$exposure_beta)

# save dictionary of id to exposure/mediators
mediator_list <- harmonised$expname$exposure
names(mediator_list) <- harmonised$expname$id.exposure
# reorder to match harmonised$exposure_beta order
mediator_list <- mediator_list[match(colnames(harmonised$exposure_beta), names(mediator_list))]

# list of exposure/mediator order in S_YY
S_YY_order <- unique(unlist(strsplit(colnames(S_XY[, c(-1, -2)]), "trait|_"))[c(FALSE, TRUE, FALSE)])

# format data for MVMR library
mvmr_data <- format_mvmr(BXGs = harmonised$exposure_beta,
    BYG = harmonised$outcome_beta,
    seBXGs = harmonised$exposure_se,
    seBYG = harmonised$outcome_se,
    RSID = rownames(harmonised$exposure_beta))

# rearrange S_YY to match order of exposure/mediators in mvmr_data
S_YY <- S_YY[match(mediator_list, S_YY_order), match(mediator_list, S_YY_order)]

# run MVMR
# robust MVMR with Q-statistic minimisation
#res1 <- qhet_mvmr(mvmr_data, S_YY, CI = TRUE, iterations = 100)
res1 <- qhet_mvmr(mvmr_data, S_YY, CI = FALSE)
saveRDS(res1, paste0(outdir, "robust_MVMR.rds"))

# MVMR with inverse variance weighted meta-analysis
res2 <- ivw_mvmr(mvmr_data)
saveRDS(res2, paste0(outdir, "ivw_MVMR.rds"))

# From https://bcgov.github.io/elucidate/reference/mean_ci.html, if you get error:
# "estimated adjustment 'a' is NA" when ci_type is set to "bca", then try again
# with more replications.