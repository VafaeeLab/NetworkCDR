#
# Takes the base graph, and removes impossible interactions
#

# We use variables and functions defined in this script
source("analysisScripts/filterNodes.r")
library("dplyr")
library("igraph")
library("pbapply")

### Credit to Jonno Bourne for union function (Stack Overflow)
# The union of two or more graphs are created. 
# The graphs may have identical or overlapping vertex sets.
union2 <- function(g1, g2, ...) {
    #Internal function that cleans the names of a given attribute
    CleanNames <- function(g, target) {
        # Get target names, find names with "_x" and remove
        gNames <- parse(text = (paste0(target, "_attr_names(g)"))) %>% eval 
        AttrNeedsCleaning <- grepl("(_\\d)$", gNames)
        StemName <- gsub("(_\\d)$", "", gNames)

        #replace attribute name for all attributes (THIS IS NOT NEEDED)
        # NewnNames <- unique(StemName[AttrNeedsCleaning])
        # for (i in NewnNames) {
        #     attr1 <- parse(text = (paste0(target, "_attr(g,'", paste0(i, "_1"), "')"))) %>% eval
        #     attr2 <- parse(text = (paste0(target, "_attr(g,'", paste0(i, "_2"), "')"))) %>% eval

        #     g <- parse(text = (paste0("set_", target, "_attr(g, i, value = ifelse(is.na(attr1), attr2, attr1))"))) %>%
        #                 eval

        #     g <- parse(text = (paste0("delete_", target, "_attr(g,'", paste0(i, "_1"), "')"))) %>% eval
        #     g <- parse(text = (paste0("delete_", target, "_attr(g,'", paste0(i, "_2"), "')"))) %>% eval
        # }
        return(g)
    }

    g <- igraph::union(g1, g2, ...) 
    for (i in c("graph", "edge", "vertex")) {
        g <- CleanNames(g, i)
    }

    return(g)
}

################################################################################

# Remove human PPI that are physically separated by subcellular location
e.filter.subcell <- function(graph) {
    # Create each subgraph and combine them
    subcell.locations <- unique(unlist(subcell$All.location))

    newg <- do.call(union2,
        pblapply(subcell.locations, function(l) {
            v.filter.subcell(graph, l)
        })
    )

    edgeDelete <- E(graph)[
        ends(graph, E(graph))[, 1] %in% V(graph)[V(graph)$type == "human-protein"]$name &
        ends(graph, E(graph))[, 2] %in% V(graph)[V(graph)$type == "human-protein"]$name
    ]
    newg <- union2(newg, delete_edges(graph, edgeDelete))
    return(newg)
}

# Trying to make more efficient version of the function
# e.filter.subcell <- function(graph) {
#     # Create each subgraph and combine them
#     subcell.locations <- unique(unlist(subcell$All.location))

#     genesInSubcell <- lapply(subcell.locations, function(l) {
#         filteredGenes <- subcell[sapply(
#         subcell$All.location, function(x) l %in% x), ]$Gene

#         filteredGenes <- filteredGenes[filteredGenes %in% V(graph)$name]
#         filteredGenes <- which(V(graph)$name %in% filteredGenes)
#         return(filteredGenes)
#     })
#     genesInSubcell <- as.data.frame(genesInSubcell)


#     newg <- graph
#     pblapply(E(graph), function(loc) {
#         lapply(subcell.locations, function(e) {
            
#             if ()
#         })
#     })
# }

# Only keep human PPI if both proteins are expressed in a tissue above a threshold
# Requires an igraph object and a dataframe with three columns: 
# Tissue, minimum nTPM and maximum nTPM
e.filter.tissue <- function(graph, ntpmCutoff = NULL) {
    if (is.null(ntpmCutoff)) {
        ntpmCutoff <- data.frame(
            tissue = unique(tissue$Tissue),
            min = rep(1, times = length(unique(tissue$Tissue))),
            max = rep(1000000, times = length(unique(tissue$Tissue)))
        )
    }

    newg <- do.call(union2, 
        pbapply(ntpmCutoff, 1, function(x) {
            v.filter.tissue(graph, x[1], x[2], x[3])
        })
    )
    
    edgeDelete <- E(graph)[
        ends(graph, E(graph))[, 1] %in% V(graph)[V(graph)$type == "human-protein"]$name &
        ends(graph, E(graph))[, 2] %in% V(graph)[V(graph)$type == "human-protein"]$name
    ]
    # newg <- union2(newg, delete_edges(graph, edgeDelete))
    return(newg)
}

# Only keep human PPI if both proteins are expressed in a cell type above a threshold
# Requires an igraph object and a dataframe with three columns: 
# Cell type, minimum nTPM and maximum nTPM
e.filter.type <- function(graph, ntpmCutoff = NULL) {
    if (is.null(ntpmCutoff)) {
        ntpmCutoff <- data.frame(
            cellType = unique(cellType$Cell.type),
            min = rep(1, times = length(unique(cellType$Cell.type))),
            max = rep(1000000, times = length(unique(cellType$Cell.type)))
        )
    }

    newg <- do.call(union2, 
        pbapply(ntpmCutoff, 1, function(x) {
            v.filter.type(graph, x[1], x[2], x[3])
        })
    )

    edgeDelete <- E(graph)[
        ends(graph, E(graph))[, 1] %in% V(graph)[V(graph)$type == "human-protein"]$name &
        ends(graph, E(graph))[, 2] %in% V(graph)[V(graph)$type == "human-protein"]$name
    ]
    newg <- union2(newg, delete_edges(graph, edgeDelete))
    return(newg)
}

# Only keep human PPI if both proteins are expressed in a cell line above a threshold
# Requires an igraph object and optionally a dataframe with three columns: 
# Cell line, minimum nTPM and maximum nTPM
e.filter.line <- function(graph, ntpmCutoff = NULL) {
    if (is.null(ntpmCutoff)) {
        ntpmCutoff <- data.frame(
            line = unique(cellLine$Cell.line),
            min = rep(1, times = length(unique(cellLine$Cell.line))),
            max = rep(1000000, times = length(unique(cellLine$Cell.line)))
        )
    }

    newg <- do.call(union2, 
        pbapply(ntpmCutoff, 1, function(x) {
            v.filter.line(graph, x[1], x[2], x[3])
        })
    )

    edgeDelete <- E(graph)[
        ends(graph, E(graph))[, 1] %in% V(graph)[V(graph)$type == "human-protein"]$name &
        ends(graph, E(graph))[, 2] %in% V(graph)[V(graph)$type == "human-protein"]$name
    ]
    newg <- union2(newg, delete_edges(graph, edgeDelete))
    return(newg)
}

################################################################################

### Sanity Checks
# g <- readRDS("clean/baseGraph.rds")
# head(V(g)) # 32902
# head(E(g)) # 673223

# newg <- e.filter.subcell(g)
# head(V(newg)) # 32902
# head(E(newg)) # 646329

# # > 1200 cell lines; Very slow
# newg <- e.filter.line(g)
# head(V(newg))
# head(E(newg))

# newg <- e.filter.type(g)
# head(V(newg))
# head(E(newg))

# newg <- e.filter.tissue(g)
# head(V(newg))
# head(E(newg))


