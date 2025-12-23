#
# Filter functions for ppi and graph
# An nTPM cutoff of 1 is used to determine if protein is expressed (https://v22.proteinatlas.org/about/assays+annotation)
#

library("jsonlite")
library("igraph")
library("dplyr")

if (!exists("filterDfLoaded")) {
    filterDfLoaded <- TRUE

    tissue <- fromJSON("clean/tissue.json")
    cellLine <- read.delim("clean/cellLine.tsv")
    cellType <- fromJSON("clean/cellType.json")
    subcell <- fromJSON("clean/subcell.json")
}

# Filter a given ppi dataframe based on an expression dataframe containing a gene,
# a type (what tissue type or cell type was measured in), and nTPM (expression level)
# For subcellular data, use `v.filter.subcell`
#
# This filters the ppi so that only proteins present in the cell type are left,
# limiting vertices to only be from the given type
v.filter <- function(ppiFull, dfExpression, type, minNTPM = 1, maxNTPM = 1000) {
    filtered <- dfExpression[dfExpression[, 2] == type, ]
    filtered <- filtered[filtered$nTPM >= minNTPM, ]
    filtered <- filtered[filtered$nTPM <= maxNTPM, ]
    filtered <- filtered[, 1]

    ppi <- ppiFull[ppiFull$P1 %in% filtered & ppiFull$P2 %in% filtered, ]
    return(ppi)
}

# Special function for subcellular data, which is not given in nTPM, but instead
# whether the protein is found or not in that area
#
# This filters the ppi so that only proteins present in the subcellular location are left,
# limiting vertices to only be from the given type
v.filter.subcell <- function(ppiFull, dfExpression = subcell, location) {
    filtered <- dfExpression[sapply(dfExpression[, 2], function(x) location %in% x), ]
    filtered <- filtered[, 1]

    ppi <- ppiFull[ppiFull$P1 %in% filtered & ppiFull$P2 %in% filtered, ]
    return(ppi)
}

# Filter a given ppi dataframe based on an expression dataframe containing a gene,
# a type (what tissue type or cell type was measured in), and nTPM (expression level)
# For subcellular data, use `e.filter.subcell`
#
# This filters edges so that only edges where both proteins are in the same cell
# type remain
e.filter <- function(ppiFull, dfExpression, minNTPM = 1, maxNTPM = 1000) {
    types <- unique(dfExpression[, 2])
    ppiRemove <- ppiFull
    for (type in types) {
        ppiKeep <- v.filter(ppiRemove, dfExpression, type, minNTPM, maxNTPM)
        ppiRemove <- anti_join(ppiRemove, ppiKeep, by = c("P1", "P2"))
    }

    ppi <- anti_join(ppiFull, ppiRemove, by = c("P1", "P2"))
    return(ppi)
}

# Special function for subcellular data, which is not given in nTPM, but instead
# whether the protein is found or not in that area
#
# This filters edges so that only edges where both proteins are in the same 
# subcellular location remain
e.filter.subcell <- function(ppiFull, dfExpression) {
    types <- unique(unlist(dfExpression[, 2]))
    ppiRemove <- ppiFull
    for (type in types) {
        ppiKeep <- v.filter.subcell(ppiRemove, dfExpression, type)
        ppiRemove <- anti_join(ppiRemove, ppiKeep, by = c("P1", "P2"))
    }

    ppi <- anti_join(ppiFull, ppiRemove, by = c("P1", "P2"))
    return(ppi)
}