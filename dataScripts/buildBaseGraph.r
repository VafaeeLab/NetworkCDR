#
# Takes the cleaned data files and builds the base graph for analysis
# 

library("jsonlite")
library("igraph")

drugs <- fromJSON("clean/drugs.json")
ppi <- fromJSON("clean/ppi.json")
disease <- fromJSON("clean/disease.json")

# Build the graph using human ppi as edges
g <- graph_from_data_frame(ppi[1:2], directed = FALSE)
V(g)$type <- "human-protein"
g <- as_directed(g, mode = "mutual")

addAdjList <- function(graph, data, node = 1, adj = 2, typeName) {
    graph <- graph + vertices(data[[node]], type = typeName)
    temp <- setdiff(unique(unlist(data[[adj]])), V(graph)$name)
    graph <- graph + vertices(temp, type = "human-protein")

    edges <- do.call(rbind, apply(data, 1, function(x) {
        cbind(x[[node]], x[[adj]])
    }))
    graph <- graph + edges(as.vector(t(edges)))
    return(graph)
}

g <- addAdjList(g, drugs, 1, 3, "drug")
g <- addAdjList(g, disease, 1, 2, "disease")

################################################################################
# Sanity Check by plotting things
type_colors <- c("human-protein" = "lightblue", 
    "drug" = "green", 
    "disease" = "red")

# Plot the subgraph with colored vertices
subgraph <- induced_subgraph(g, 
        unlist(ego(g, order = 1, nodes = "CHEMBL1743017")))
plot(
    subgraph,
    vertex.color = type_colors[V(subgraph)$type],
    main = "Subgraph with Vertices Colored by 'Type'"
)

subgraph <- induced_subgraph(g, 
    unlist(ego(g, order = 1, 
        nodes = "46, XX Testicular Disorders of Sex Development")))
plot(
    subgraph,
    vertex.color = type_colors[V(subgraph)$type],
    main = "Subgraph with Vertices Colored by 'Type'"
)

subgraph <- induced_subgraph(g, 
    unlist(ego(g, order = 1, 
        nodes = "CYP51A1")))
plot(
    subgraph,
    vertex.color = type_colors[V(subgraph)$type],
    main = "Subgraph with Vertices Colored by 'Type'"
)
################################################################################

saveRDS(g, "clean/baseGraph.rds")
