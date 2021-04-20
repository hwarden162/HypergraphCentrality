library(igraph)
library(reticulate)

use_virtualenv("default")
source_python("/Users/hugh/Documents/University/Maths/Year4/Project/RCode/CentralityCode/HypergraphBC.py")

load("/Users/hugh/Documents/University/Maths/Year4/Project/RCode/CentralityCode/Data/Constructions.RData")

hypergaph.dual_adj_mat <- get_hyperedge_dual_adj_mat(ppi.inc_mat)

hypergaph.dual_adj_mat <- matrix(as.numeric(hypergaph.dual_adj_mat > 0), ncol = ncol(hypergaph.dual_adj_mat))
diag(hypergaph.dual_adj_mat) <- 0

hypergraph.bc <- get_hypergraph_betweenness(hypergaph.dual_adj_mat, ppi.inc_mat)
names(hypergraph.bc) <- c("all", "nt", "act")

graph.bc <- betweenness(ppi.graph)

bc.data <- data.frame(protein = rownames(ppi.inc_mat), gdeg = degree(ppi.graph), ncomp = apply(ppi.inc_mat, 1, sum),graph = graph.bc, hall = hypergraph.bc$all, hnt = hypergraph.bc$nt, hact = hypergraph.bc$act)

save(bc.data, file = "/Users/hugh/Documents/University/Maths/Year4/Project/RCode/CentralityCode/Data/Betweenness.RData")
