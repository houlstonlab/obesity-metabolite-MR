###
# Estimate covariance matrices using metaCCA and perform sensitivity tests for MVMR.
###

library(metaCCA)
library(tidyverse)
library(MVMR)

outdir <- "/path/to/your/output/directory/"

# estimate S_YY using metaCCA
S_XY <- read.table(paste0(outdir, "collated_harmonised.tsv"),
                    header = TRUE, row.names = 1, stringsAsFactors = TRUE)
S_YY <- estimateSyy(S_XY = S_XY)

# load harmonised MVMR data
harmonised <- readRDS(paste0(outdir, "cancer_harmonised_corrected.tsv"))

# save order of obesity measure/mediators
x_list <- colnames(harmonised$exposure_beta)

# save dictionary of id to exposure/mediators
mediator_list <- harmonised$expname$exposure
names(mediator_list) <- harmonised$expname$id.exposure
# reorder to match harmonised$exposure_beta order
mediator_list <- mediator_list[match(colnames(harmonised$exposure_beta), names(mediator_list))]

# list of exposure/mediator order in S_XY and S_YY
S_XY_order <- unlist(strsplit(colnames(S_XY[, c(-1, -2)]), "trait"))[c(FALSE, TRUE)]
S_YY_order <- unique(unlist(strsplit(colnames(S_XY[, c(-1, -2)]), "trait|_"))[c(FALSE, TRUE, FALSE)])

# format data for MVMR library
mvmr_data <- format_mvmr(BXGs = harmonised$exposure_beta,
    BYG = harmonised$outcome_beta,
    seBXGs = harmonised$exposure_se,
    seBYG = harmonised$outcome_se,
    RSID = rownames(harmonised$exposure_beta))

# rearrange S_XY and S_YY to match order of exposure/mediators in mvmr_data
S_XY <- S_XY[, c(c(1, 2), 2 + unlist(lapply(mediator_list, function(x) grep(x, S_XY_order))))]
S_YY <- S_YY[match(mediator_list, S_YY_order), match(mediator_list, S_YY_order)]

# keep only SNPs that are in mvmr_data
S_XY <- S_XY %>% filter(rownames(S_XY) %in% rownames(mvmr_data))
# and reorder
S_XY <- S_XY[match(rownames(mvmr_data), rownames(S_XY)), ]

# estimate covariance matrices from phenotypic correlation matrix
# (second argument is a table of standard errors)
cov_mat <- phenocov_mvmr(S_YY, S_XY[, seq(4, length(S_XY), 2)])

# test for weak instruments
sres <- strength_mvmr(mvmr_data, gencov = cov_mat)
saveRDS(sres, paste0(outdir, "sres.rds"))

# test for horizontal pleiotropy
pres <- pleiotropy_mvmr(mvmr_data, gencov = cov_mat)
saveRDS(pres, paste0(outdir, "pres.rds"))