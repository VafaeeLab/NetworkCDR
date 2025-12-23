#
# Clean the drug combinations data and convert the IDs
#

library(dplyr)
library(tidyr)
library(jsonlite)

# Get maps to ensure that drug combinations are given in usable IDs
mapID <- read.csv("rawdata/drug-mappings.tsv", header = TRUE, sep = "\t")
mapID <- mapID[c("drugbankId", "chembl_id")] %>% filter(chembl_id != "null")

combs <- read.csv("CDCDB/all_combs_unormalized.csv", header = TRUE)
conds <- read.csv("CDCDB/conditions_df.csv", header = TRUE)

# Change column name and then remove duplicates and NAs
names(combs)[names(combs) == "source_id"] <- "nct_id"
combs <- combs[combs$nct_id != "N/A", ]
combs <- distinct(combs, nct_id, .keep_all = TRUE)

knownCombs <- inner_join(combs, conds, by = "nct_id")
knownCombs <- knownCombs[c("nct_id", "condition", "drugbank_identifiers")]

# Get drugbank and chembl IDs
knownCombs$drug_ids <- lapply(knownCombs$drugbank_identifiers, function(x) {
    x <- gsub("\\[|\\]", "", x)
    x <- gsub("\"", "", x)
    strsplit(x, ", ")[[1]]
})

knownCombs <- knownCombs[sapply(knownCombs$drug_ids, function(x) !("-1" %in% x)), ]
knownCombs$chembl <- lapply(knownCombs$drug_ids, function(x) {
    lapply(x, function(y) {
        ifelse(y %in% mapID$drugbankId, mapID[mapID$drugbankId == y, "chembl_id"], "")
    })
})
knownCombs <- knownCombs[sapply(knownCombs$chembl, function(x) !("" %in% x)), ]
# tail(knownCombs) # Sanity check

saveRDS(knownCombs, file = "clean/verifiedCombinations.RData")

