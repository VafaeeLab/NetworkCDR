#
# Get useful statistics and explore data for useful information
# 

library("igraph")
library("ggplot2")

g <- readRDS("clean/baseGraph.rds")
length(V(g)[[type == "drug"]])
length(V(g)[[type == "human-protein"]])
length(V(g)[[type == "disease"]])

# How does the degree distribution change with filtering
df <- as.data.frame(as.numeric(degree(g)))
colnames(df) <- c("degree")
ggplot(data = df, aes(x = degree)) +
    geom_histogram(color = "black", 
        fill = "white", 
        binwidth = 5) +
    xlim(0, 100)

gFil <- e.filter.subcell(gFil)
gFil <- e.filter.tissue(gFil)
gFil <- e.filter.type(gFil)
degree_distribution(gFil)