###
# Find associations for SNPs in the GWAS Catalog.
###

library(httr)
library(jsonlite)
library(tidyverse)

iv_tibble_file <- "/path/to/metabolite/exposures/iv_list.txt"

iv_tibble <- read_delim(iv_tibble_file)
iv_tibble$gwas_catalog <- NA
print(paste("Number of SNPs to query:", nrow(iv_tibble)))

get_assocs <- function(rsid) {
    # get associations for a given rsid
    res <- GET(paste0("https://www.ebi.ac.uk/gwas/rest/api/singleNucleotidePolymorphisms/", rsid, "/associations"))
    Sys.sleep(0.2) # sleep to avoid rate limiting

    # if not OK, skip
    if (res$status_code != 200) return(NULL)

    # clean data returned from API
    data <- fromJSON(rawToChar(res$content))
    names(data) <- gsub("_", "", names(data))
    names(data$embedded$associations) <- gsub("_", "", names(data$embedded$associations))

    # get list of traits associated with SNP
    trait_links <- data$embedded$associations$links$efoTraits

    trait_string <- list()
    for (i in seq_along(trait_links[, ])) {
        # get trait data
        res <- GET(trait_links[i, ])
        Sys.sleep(0.2) # sleep to avoid rate limiting
        if (res$status_code != 200) next
        trait_data <- fromJSON(rawToChar(res$content))
        names(trait_data) <- gsub("_", "", names(trait_data))

        # since different studies put the trait name in different fields,
        # we loop through the possible fields looking for a non-empty trait name
        if (all(is.na(data$embedded$associations[i, ]$pvalueDescription))) {
            data$embedded$associations[i, ]$pvalueDescription <- ""
            if (any("metabolite measurement" %in% trait_data$embedded$efoTraits$trait)) {
                res <- GET(data$embedded$associations$links$study[i, ])
                Sys.sleep(0.2) # sleep to avoid rate limiting
                if (res$status_code != 200) next
                study_data <- fromJSON(rawToChar(res$content))
                trait_data$embedded$efoTraits$trait <- paste(trait_data$embedded$efoTraits$trait, study_data$diseaseTrait$trait, sep = " ")
            } else if (any("blood metabolite measurement" %in% trait_data$embedded$efoTraits$trait)) {
                res <- GET(data$embedded$associations$links$study[i, ])
                Sys.sleep(0.2) # sleep to avoid rate limiting
                if (res$status_code != 200) next
                study_data <- fromJSON(rawToChar(res$content))
                trait_data$embedded$efoTraits$trait <- paste(trait_data$embedded$efoTraits$trait, study_data$diseaseTrait$trait, sep = " ")
            }
        }

        if (all(is.na(trait_data$embedded$efoTraits$trait))) {
            trait_data$embedded$efoTraits$trait <- ""
        }

        # append trait name and p-value description
        trait_string <- c(trait_string, paste(paste(trait_data$embedded$efoTraits$trait, sep = " "),
                        paste(data$embedded$associations[i, ]$pvalueDescription, sep = " "), sep = " "))
    }
    if (length(trait_string) == 0) return(NA)
    return(paste(trait_string, collapse = " + "))
}

for (i in seq_len(nrow(iv_tibble))) {
    if (i %% 100 == 0) print(paste("SNP number:", i))
    assocs <- get_assocs(iv_tibble$rsid[i])
    if (is.null(assocs)) next
    iv_tibble$gwas_catalog[i] <- assocs
}

iv_tibble_file <- gsub(".txt", "", iv_tibble_file)
write_tsv(iv_tibble, paste0(iv_tibble_file, "_gwas_catalog.txt"))
