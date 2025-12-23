#
# Testing the effectiveness of the network-based repositioning method using known
# and verified pairwise drug combinations
#

library("jsonlite")
library("igraph")
library("tidyr")
library("dplyr")
library("Metrics")
library("pbapply")

source("analysisScripts/buildGraph.r")
source("analysisScripts/filters.r")
source("analysisScripts/findDistance.r")
source("analysisScripts/findDrugs.r")
source("analysisScripts/visualisation.r")


g <- readRDS("clean/baseGraph.rds")
knownCombs <- readRDS("clean/verifiedCombinations.RData")

sum(tolower(gsub("[[:punct:] ]", "", knownCombs$condition)) %in% tolower(gsub("[[:punct:] ]", "", V(g)$name)))
length(knownCombs$condition)

# Need to get this working
library(stringdist)
D <- stringdistmatrix(knownCombs$condition, V(g)$name, method = "jw")
bestMatch <- apply(D, 1, which.min)

df <- data.frame(
  knownCombCond = knownCombs$condition,
  bestGraphMatch = V(g)$name[bestMatch],
  distance = D[cbind(seq_along(knownCombs$condition), bestMatch)]
)


# Verify this works first
checkDrugComb(g, knownCombs$condition[[2]], knownCombs$chembl[[2]])
knownCombs$chembl[[2]]
V(g)["CHEMBL385517"] # Need to associate alternate forms of drugs together 

knownCombs$condition[[1]]
V(g)["Trigger Finger"] # Need to fuzzy search or use some other search method

knownCombs$myResults %>% 
    lapply(function(x) {
        if (is.data.frame(x)) return(FALSE)
        return(x == "Error: No disease")
        return(TRUE)
    }) %>%
    unlist() %>% 
    sum()


################################################################################

# OLD METHOD (KEPT AS REFERENCE)

checkDrugDist <- function(graph, drugPairs, distFunc) {
    unlist(apply(drugPairs, MARGIN = 1, function(comb) {
        nodeset1 <- neighbors(graph, comb[[1]], mode = "out")
        nodeset2 <- neighbors(graph, comb[[2]], mode = "out")
        return(distFunc(graph, nodeset1, nodeset2))
    }))
}

# Since the drug combinations do not say what disease they are for, 
# we will find the closest diseases
g <- readRDS("clean/baseGraph.rds")
knownCombs$diseases <- apply(knownCombs, MARGIN = 1, function(drugs) {
    if (!(drugs[[1]] %in% V(g)$name && drugs[[2]] %in% V(g)$name)) {
        return(NULL)
    }
    dis <- ego(g, order = 2, nodes = c(drugs[[1]], drugs[[2]]), mode = "all")
    dis <- dis[[1]][dis[[1]] %in% dis[[2]]]
    dis <- dis[dis$type == "disease"]
    return(dis$name)
})

lengths <- sapply(knownCombs$diseases, function(x) length(x))
knownCombs <- knownCombs[lengths != 0, ]

##### Find predicted drug pairs that display complementary exposure
knownCombs$topDist <- checkDrugDist(g, knownCombs, topDist)
knownCombs$overlapDist <- checkDrugDist(g, knownCombs, overlapDist)

nrow(knownCombs) # 595 drug combinations
knownCombs %>% filter(overlapDist >= 0) %>% nrow() # 505 drug combinations
knownCombs %>% filter(topDist >= 0) %>% nrow() # 352 drug combinations


##### We can now repeat but with filtered edges
gFil <- readRDS("clean/baseGraph.rds")
gFil <- e.filter.subcell(gFil)
gFil <- e.filter.tissue(gFil)
gFil <- e.filter.type(gFil)

drugs <- fromJSON("clean/drugs.json")
ppi <- fromJSON("clean/ppi.json")
disease <- fromJSON("clean/disease.json")

ppiFil <- e.filter.subcell(ppi, subcell) %>%
    e.filter(tissue) %>%
    e.filter(cellType)
gFil <- buildGraph(ppiFil, drugs, disease)

filtered <- knownCombs
filtered$topDist <- checkDrugDist(gFil, filtered, topDist)
filtered$overlapDist <- checkDrugDist(gFil, filtered, overlapDist)

nrow(filtered) # x drug combinations
filtered %>% filter(overlapDist >= 0) %>% nrow() # x drug combinations
filtered %>% filter(topDist >= 0) %>% nrow() # x drug combinations

##### We can check for false positives using random drug pairs
# Utilising method from https://www.nature.com/articles/s41467-019-09186-x

findAUC <- function(graph, drugPairs) {
    genPairs <- t(replicate(nrow(drugPairs), 
        sample(V(graph)[V(graph)$type == "drug"]$name, size = 2, replace = FALSE)))
    
    genPairsTopDist <- checkDrugDist(graph, drugPairs, topDist)

    finalAUC <- auc(c(rep(TRUE, nrow(drugPairs)), rep(FALSE, nrow(drugPairs))), 
        c(drugPairs$topDist >= 0, genPairsTopDist < 0))
    
    return(finalAUC)
}

noFiltMean <- mean(pbreplicate(100, findAUC(g, knownCombs))) # 0.5915966
filtMean <- mean(pbreplicate(100, findAUC(gFil, filtered)))   


################################################################################
################################################################################
# After acquiring known combinations from clinical trials site, I can test my methods

source("analysisScripts/filterNodes.r")
source("analysisScripts/filterEdges.r")
source("analysisScripts/findDistance.r")
source("analysisScripts/findDrugs.r")

g <- readRDS("clean/baseGraph.rds")

# Get combinations that have been verified
knownCombs <- readRDS("knownCombs.rds")
knownCombs$myResults <- apply(knownCombs, MARGIN = 1, function(x) {
    checkDrugComb(g, x$condition, x$chembl)
})


##### We can now repeat but with filtered edges
gFil <- readRDS("clean/baseGraph.rds")
gFil <- e.filter.subcell(gFil)
gFil <- e.filter.tissue(gFil)
gFil <- e.filter.type(gFil)

filtered <- readRDS("knownCombs.rds")
filtered$myResults <- apply(filtered, MARGIN = 1, function(x) {
    checkDrugComb(g, x$condition, x$chembl)
})




