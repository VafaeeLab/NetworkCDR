#
# Find topological distances between network nodes
#

library("igraph")
library("dplyr")
library("pbapply")
library("parallel")

if (!exists("GOSemSimLoaded")) {
    GOSemSimLoaded <- TRUE

    library("GOSemSim")
    library("org.Hs.eg.db")
    goDataCC <- godata(annoDb = "org.Hs.eg.db", ont="CC", computeIC=FALSE)
    goDataMF <- godata(annoDb = "org.Hs.eg.db", ont="MF", computeIC=FALSE)
    goDataBP <- godata(annoDb = "org.Hs.eg.db", ont="BP", computeIC=FALSE)
}

# Find the average shortest distances between two sets of nodes
# Implements formula: s_AB = d_AB - (d_AA + d_BB)/2
# Where d_XY is the sum of shortest distances from nodes in X to any node in Y 
# and nodes in Y to any node in X, divided by the number of nodes
topDist <- function(graph, nodeset1, nodeset2) {
    changeInfTo <- length(V(graph)) # May change if better method found
    
    d <- distances(graph, v = nodeset1, to = nodeset2, mode = "out")
    d[d == Inf] <- changeInfTo 
    d_AB <- sum(apply(d, 1, min)) + sum(apply(d, 2, min))
    d_AB <- d_AB / (nrow(d) + ncol(d))

    d <- distances(graph, v = nodeset1, to = nodeset1, mode = "out")
    diag(d) <- Inf # Remove the zeroes in diagonal
    d[d == Inf] <- changeInfTo
    d_AA <- sum(apply(d, 1, min))
    d_AA <- d_AA / nrow(d)

    d <- distances(graph, v = nodeset2, to = nodeset2, mode = "out")
    diag(d) <- Inf
    d[d == Inf] <- changeInfTo
    d_BB <- sum(apply(d, 1, min))
    d_BB <- d_BB / nrow(d)

    return(d_AB - (d_AA + d_BB)/2)
}

# Alternative distance function given by `shortest separation - overlap`
overlapDist <- function(graph, nodeset1, nodeset2) {
    changeInfTo <- length(V(graph))

    d <- min(distances(graph, v = nodeset1, to = nodeset2, mode = "out"))
    ans <- ifelse(d == Inf, changeInfTo, d - sum(nodeset1 %in% nodeset2))
    return(ans)
}

# Function to find functional proximity between two sets of nodes
# Ensure go data is preprocessed for function to work (should occur when you run this script)
# `out` can be "all", "max", "min", or "mean"; `out` is ignored if `type` is not "all"
# `type` can be "CC" (Cellular Component), "MF" (Molecular Function), "BP" (Biological Process), or "all"
funcDist <- function(graph, nodeset1, nodeset2, type = "BP", out = "min") {
    set1 <- select(org.Hs.eg.db,
        keys = names(nodeset1),
        columns = "ENTREZID",
        keytype = "SYMBOL")[, "ENTREZID"] %>% suppressMessages()
    set2 <- select(org.Hs.eg.db,
        keys = names(nodeset2),
        columns = "ENTREZID",
        keytype = "SYMBOL")[, "ENTREZID"] %>% suppressMessages()

    if (type == "all") {
        simCC <- 1 - mgeneSim(c(set1, set2), goDataCC, verbose = FALSE)[set1, set2]
        simMF <- 1 - mgeneSim(c(set1, set2), goDataMF, verbose = FALSE)[set1, set2]
        simBP <- 1 - mgeneSim(c(set1, set2), goDataBP, verbose = FALSE)[set1, set2]
    } else if (type == "CC") {
        simCC <- 1 - mgeneSim(c(set1, set2), goDataCC, verbose = FALSE)[set1, set2]
        return(mean(simCC))
    } else if (type == "MF") {
        simMF <- 1 - mgeneSim(c(set1, set2), goDataMF, verbose = FALSE)[set1, set2]
        return(mean(simMF))
    } else if (type == "BP") {
        simBP <- 1 - mgeneSim(c(set1, set2), goDataBP, verbose = FALSE)[set1, set2]
        return(mean(simBP))
    } else {
        print(paste("Invalid value for type: ", type))
        return()
    }

    res <- c(mean(simCC), mean(simMF), mean(simBP))
    if (out == "all") {
        return(res)
    } else if (out == "mean") {
        return(mean(res))
    } else if (out == "max") {
        return(max(res))
    } else if (out == "min") {
        return(min(res))
    } else {
        print(paste("Invalid value for out: ", out))
        return()
    }
}

################################################################################

# Given a graph and disease name, finds nearby drugs and returns their distance values
# `minSep` and `maxSep` are the minimum and maximum separation, 
# representing the number of nodes that can appear between the disease and drugs
# `diseaseDistShift` is the amount to subtract from drug-disease distances. 
# This would allow further drugs to be included, which can be finetuned to include 
# drugs that affect proteins in the neighbourhood of the disease.
findDrugDist <- function(graph, disease, minSep = 0, maxSep = 0, distFun = topDist, diseaseDistShift = maxSep) {
    v <- ego(graph, order = maxSep+2, nodes = disease, mode = "all")[[1]]
    v <- v %m% ego(graph, order = minSep+1, nodes = disease, mode = "all")[[1]]

    drugs <- v[v$type == "drug"]
    drug_names <- names(drugs)
    cat(paste("Found", length(drug_names), "drugs", "\n"))

    diseaseCluster <- neighbors(graph, disease, mode = "out")
    drugCluster <- lapply(drugs, function(x) neighbors(graph, x, mode = "out"))

    cat("Calculating disease-drug distance\n")
    diseaseDist <- unlist(pblapply(drug_names,
        function(x) distFun(graph, diseaseCluster, drugCluster[[x]])))

    diseaseDist <- diseaseDist - diseaseDistShift
    drug_names <- drug_names[diseaseDist <= 0]

    cat(paste("Found", length(drug_names), "useful drugs\n"))
    if (length(drug_names) == 0) {
        return(diseaseDist)
    }

    diseaseDist <- diseaseDist[diseaseDist <= 0]

    cat("Calculating drug-drug distance\n")
    distM <- do.call(rbind, pblapply(drug_names, function(i) {
        mclapply(drug_names, 
            function(j) distFun(graph, drugCluster[[i]], drugCluster[[j]])
        )
    }))
    dimnames(distM) <- list(drug_names, drug_names)

    return(data.frame(cbind(diseaseDist, distM)))
}

# Given a graph, disease, and drug combination list provides distances based on 
# the distance function `distFun` 
checkDrugComb <- function(graph, disease, drugComb, distFun = topDist) {
    drugComb <- unlist(drugComb)

    if (!(disease %in% V(graph)$name)) {
        print("Error: No disease")
        return("Error: No disease")
    }

    if (lapply(drugComb, function(x) !(x %in% V(graph)$name)) %>%
        unlist() %>%
        any()) {
            print("Error: No drug")
        return("Error: No drug")
    }

    diseaseCluster <- neighbors(graph, disease, mode = "out")
    drugCluster <- lapply(V(graph)[drugComb], function(x) neighbors(graph, x, mode = "out"))

    cat("Calculating disease-drug distance\n")
    diseaseDist <- unlist(pblapply(drugComb,
        function(x) distFun(graph, diseaseCluster, drugCluster[[x]])))

    cat("Calculating drug-drug distance\n")
    distM <- do.call(rbind, pblapply(drugComb, function(i) {
        mclapply(drugComb, 
            function(j) distFun(graph, drugCluster[[i]], drugCluster[[j]])
        )
    }))
    dimnames(distM) <- list(drugComb, drugComb)

    return(data.frame(cbind(diseaseDist, distM)))
}

### Sanity Check
# g <- readRDS("clean/baseGraph.rds")
# findDrugDist(g, "ACROMESOMELIC DYSPLASIA, MAROTEAUX TYPE", 0, 0)
# findDrugDist(g, "ACROMESOMELIC DYSPLASIA, MAROTEAUX TYPE", 0, 0, overlapDist)

# checkDrugComb(g, chosenDisease, drugCombs[[2]][[2]])

# s1 <- neighbors(g, "CHEMBL1089636", mode = "out")
# s2 <- neighbors(g, "CHEMBL1364551", mode = "out")
# funcDist(g, s1, s2, out = "all", type = "all")
# funcDist(g, s1, s2)


# findDrugDist(g, "ALZHEIMER DISEASE 2", 0, 0, funcDist, diseaseDistShift = 0.6)
# checkDrugComb(g, chosenDisease, drugCombs[[2]][[1]], funcDist)
