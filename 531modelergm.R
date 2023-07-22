# Set the working directory
setwd("F:/ergm")

# Install the intergraph package if it is not installed
if (!requireNamespace("intergraph", quietly = TRUE)) {
  install.packages("intergraph")
}

# Load required packages
# Load the parallel package
library(parallel)
library(network)
library(Rglpk)
library(sna)
library(igraph)
library(ergm)
library(ergm.count)
#library(ergm.wergm)
library(texreg)
library("readxl")
library("statnet")
library("intergraph")


# Read the datasets
data <- read.csv("All(2019).csv", header = TRUE)

# Convert the edge list to an igraph object
g <- graph_from_data_frame(data, directed = FALSE)

# Convert the igraph object to a network object with edge weights
net_weighted <- asNetwork(g)

# Set edge weights in the network object
network::set.edge.attribute(net_weighted, "weights", data$weight)

# Create a matrix of edge weights
edge_weights_matrix <- as.matrix(net_weighted %e% "weights")

# Weight and edges
model1 <- ergm(net_weighted ~ edgecov(edge_weights_matrix) + edges(), 
               control = control.ergm(MCMC.samplesize = 80000, 
                                      MCMC.burnin = 1000000, 
                                      MCMC.interval = 10000, 
                                      MCMLE.density.guard = 100000, 
                                      seed = 567))

summary(model1)

# Model 2:Read the 2019-1.csv dataset
node_attr <- read.csv("All(2019-11).csv")

# Create a network object with node covariates
net_nodecov <- network(net_weighted, vertex.attr = node_attr, directed = FALSE)

# Fit a W-ERGM model with nodeGDP only
model2 <- ergm(net_nodecov ~ nodecov("nodeType")+nodecov("nodeGDP")+nodecov("nodeTheme")+nodecov("nodeCode"), control = control.ergm(seed = 567))

summary(model2)


# Exogenous network effects
# Read in the data
distance_matrix <- read.csv("dist_matrix.csv", header = FALSE)

# Remove duplicated rows
distance_matrix <- unique(distance_matrix)

# Create the network object
net_binary <- network(distance_matrix, directed = FALSE)

# Convert the network object to an edgelist
edgelist1 <- as.edgelist(net_binary)

contig_matrix <- as.matrix(read.csv("contig_matrix.csv", header = FALSE))

# Set the diagonal elements to 0
diag(contig_matrix) <- 0

# Create the network object
net_binary <- network(contig_matrix, directed = FALSE)

# Convert the network object to an edgelist
edgelist2 <- as.edgelist(net_binary)

# Combine the adjacency matrices
combined_matrix <- distance_matrix + contig_matrix

# Create the combined network object
net_combined <- network(combined_matrix, directed = FALSE)

# Create the combined edgelist
edgelist_combined <- as.edgelist(net_combined)

# Fit a W-ERGM model with binary  distance
model3 <- ergm(net_binary ~ edgecov(edgelist1) + edgecov(edgelist2), control = control.ergm(seed = 567))

summary(model3)

# Fit an ERGM model with edgelist1, edgelist2, edge_weights_matrix, nodeType and nodeTheme
model14 <- ergm(net_combined ~ edgecov(edgelist1) + edgecov(edgelist2) + edgecov(edge_weights_matrix) + edges() ,
                control = control.ergm(seed = 567))
# Summary of the model
summary(model14)

# Actor Attribute Model
# Model 2:Read the 2019-1.csv dataset
node_attr <- read.csv("All(2019-11).csv")

# Create a network object with node covariates
net_nodecov <- network(net_weighted, vertex.attr = node_attr, directed = FALSE)

# Create the combined network object with node attributes
net_combined <- network(combined_matrix, directed = FALSE, vertex.attr = node_attr)
# Fit a W-ERGM model with nodeGDP only
model5 <- ergm(net_combined ~ edgecov(edge_weights_matrix) +edges() + nodecov("nodeType")+nodecov("nodeGDP")+nodecov("nodeTheme")+nodecov("nodeCode") , control = control.ergm(seed = 567))
# Summary of the model
summary(model5)


model6 <- ergm(net_combined ~ edgecov(edgelist1) + edgecov(edgelist2) + edgecov(edge_weights_matrix) +edges() + nodecov("nodeType")+nodecov("nodeGDP")+nodecov("nodeTheme")+nodecov("nodeCode") , control = control.ergm(seed = 567))
# Summary of the model
summary(model6)