#
# Takes the data files and makes their formats consistent.
# Use preferred name for all proteins and chembl for drugs
# 

library("jsonlite")
library("STRINGdb")
library("dplyr")
library("parallel")
library("pbapply")
library("tidyr")
library("stringr")

# Get protein database to make all names the same
string_db <- STRINGdb$new(
    version = "12",
    species = 9606,
    score_threshold = 900
)

proteindb <- string_db$get_aliases()

preferredName <- function(x) {
    stringid <- proteindb$STRING_id[match(x, proteindb$alias)]
    return(proteindb$alias[match(stringid, proteindb$STRING_id)])
}

################################################################################

# Fix PPI data
ppi <- read.table("rawdata/stringdb-human-known_ppi-confident.txt", 
    stringsAsFactors = FALSE, 
    col.names = c("P1", "P2", "experimental", "database"))
ppi$P1 <- preferredName(ppi$P1)
ppi$P2 <- preferredName(ppi$P2)

# Fix drug data - Note the lapply is quite slow
drugs <- fromJSON(txt="rawdata/opentargets-drug.json")
drugs <- drugs[c("id", "name", "linkedTargets")]
drugs$linkedTargets <- drugs$linkedTargets[["rows"]]
drugs$linkedTargets <- pblapply(drugs$linkedTargets, function(x) {
    return(mclapply(x, preferredName, mc.cores = detectCores()))
})
drugs <- drugs[!is.na(drugs$linkedTargets), ]

# Fix disease data
disgenet <- readLines("rawdata/DisGeNET.txt")
disgenet <- strsplit(disgenet, "\t")
disgenet <- data.frame(
    disease = sapply(disgenet, function(x) x[1]),
    genes = I(lapply(disgenet, function(x) x[-1]))
)
disgenet <- disgenet %>% 
    rowwise() %>% 
    mutate(genes = list(genes[genes != ""])) %>% 
    ungroup()

# Fix gene tissue data
tissue <- read.delim("rawdata/rna_tissue_consensus.tsv")
tissue$Gene <- tissue$Gene.name
tissue$Gene.name <- NULL

# Fix gene cell line data
cellLine <- read.delim("rawdata/rna_celline.tsv")
cellLine$Gene <- cellLine$Gene.name
cellLine <- cellLine[c("Gene", "Cell.line", "nTPM")]

# Fix gene cell type data
celltype <- read.delim("rawdata/rna_single_cell_type.tsv")
celltype$Gene <- celltype$Gene.name
celltype <- celltype[c("Gene", "Cell.type", "nTPM")]

# Fix gene cell subcellular data
subcell <- read.delim("rawdata/subcellular_location.tsv")
subcell$Gene <- subcell$Gene.name
subcell$GO.id <- str_remove_all(subcell$GO.id, " \\(GO:[0-9]*\\)")
subcell$GO.id <- str_remove_all(subcell$GO.id, " \\(\\)")
subcell$All.location <- strsplit(subcell$GO.id, ";")
subcell$Main.location <- strsplit(subcell$Main.location, ";")
subcell <- subcell[c("Gene", "All.location", "Main.location")]

################################################################################
# Write as clean data to json
write(toJSON(drugs), "clean/drugs.json")
write(toJSON(ppi), "clean/ppi.json")
write(toJSON(disgenet), "clean/disease.json")
write(toJSON(tissue), "clean/tissue.json")
write.table(cellLine, "clean/cellLine.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
write(toJSON(celltype), "clean/cellType.json")
write(toJSON(subcell), "clean/subcell.json")
