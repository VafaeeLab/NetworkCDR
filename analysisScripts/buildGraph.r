#
# Given a filtered ppi dataset, build a graph
#

library("jsonlite")
library("igraph")

buildGraph <- function(ppi, drugs, disease) {
    # Create a graph containing all protein-protein interactions as bidirectional edges
    g <- graph_from_data_frame(ppi[1:2], directed = FALSE)
    V(g)$type <- "human-protein"
    g <- as_directed(g, mode = "mutual")

    # Add nodes and edges to the graph based on an adjacency list
    # `graph` is the graph to add to; `data` is the adjacency list as a dataframe
    # `node` is the column in `data` which contains the node name
    # `adj` is the column in `data` which contains the adjacent protein
    # `typeName` is the type given to the nodes added through this function
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
}


