#!/bin/bash

###
# Script to run MVMR analysis for an obesity measure with metabolite mediators
# and a cancer outcome.
###

# select only SNPs that are IVs for either the exposure or a mediator
Rscript collate_exposures.r M00001 M00002

# reclump IVs selected above
Rscript MVMR_clump.r M00001 M00002

# harmonise full GWAS summary stats for exposure and mediators for use in metaCCA
Rscript MVMR_sensitivity_collate.r obesity_measure M00001 M00002

# remove palindromic SNPs that are ambiguous
Rscript correct_MVMR_SNPs.r cancer_name

# run MVMR sensitivity tests
Rscript MVMR_sensitivity_tests.r

# run MVMR
Rscript robust_MVMR.r