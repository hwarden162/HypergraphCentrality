library(igraph)

setwd("/Users/hugh/Documents/University/Maths/Year4/Project/RCode/CentralityCode")

mips.full <- read.csv("/Users/hugh/Documents/University/Maths/Year4/Project/RCode/CentralityCode/Data/mips.txt", sep = "\t")
mips.full <- mips.full[1:217,]
mips.edgelist <- mips.full[,1:2]

ppi.full <- as.data.frame(read.csv("/Users/hugh/Documents/University/Maths/Year4/Project/RCode/CentralityCode/Data/BIOGRID-ORGANISM-Saccharomyces_cerevisiae_S288c-4.3.195.tab3.txt", sep = "\t"))
ppi.edgelist <- ppi.full[,6:7]

mips.proteins <- unique(c(mips.edgelist[,1], mips.edgelist[,2]))
ppi.proteins <- unique(c(ppi.edgelist[,1], ppi.edgelist[,2]))
ppi.edgelist <- unique(ppi.edgelist[which((ppi.edgelist[,1] %in% mips.proteins)&(ppi.edgelist[,2] %in% mips.proteins)),])

ppi.graph <- graph_from_edgelist(as.matrix(ppi.edgelist), directed = FALSE)
ppi.adj_mat <- as.matrix(as_adjacency_matrix(ppi.graph))
rnames <- rownames(ppi.adj_mat)
ppi.adj_mat <- matrix(as.numeric(ppi.adj_mat > 0), ncol = ncol(ppi.adj_mat))
diag(ppi.adj_mat) <- 0
ppi.graph <- graph_from_adjacency_matrix(ppi.adj_mat, mode = "undirected", diag = FALSE)

ppi.cliques <- max_cliques(ppi.graph, min = 2)

ppi.inc_mat <- matrix(0, nrow = length(rnames), ncol = length(ppi.cliques))
rownames(ppi.inc_mat) <- rnames

for (i in 1:length(ppi.cliques)){
  ppi.inc_mat[ppi.cliques[[i]],i] <- 1
}


save(ppi.graph, ppi.inc_mat, file = "/Users/hugh/Documents/University/Maths/Year4/Project/RCode/CentralityCode/Data/Constructions.RData")


