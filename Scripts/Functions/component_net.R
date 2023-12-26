# --------------------- component structure --------------------

# Author: Kate P Maia
# Date: 06/2020 (Checked: 12/2023)

# component_str: creates dataframe with information on the network-level component structure. 
# input: M (incidence interaction matrix [I, J])
# output: df with number of components (NComp), size of the largest component (CSize) and number of components with the largest size (nLarge).

library(igraph)

# --------------------------------------------------------------

component_str <- function(M) {
  
  # Adjacency from incidence matrix
  nA <- nrow(M); nP <- ncol(M)
  q2 <- q4 <- M 
  q1 <- mat.or.vec(nA, nA); q3 <- mat.or.vec(nP, nP)
  upper <- cbind(q1, q2); lower <- cbind(t(q4), q3)
  colnames(upper) <- colnames(lower) <- c(rownames(M), colnames(M))
  A <- as.matrix(rbind(upper, lower))
  rownames(A) <- c(rownames(M), colnames(M))
  if (any(rownames(A) != colnames(A))) {print("ERROR IN COMPONENT_STR: 1")}
  
  g <- graph_from_adjacency_matrix(A, "undirected")
  clu <- igraph::components(g)
  
  NComp <- clu$no # number of components
  CSize <- max(clu$csize) # size of largest component
  largclu <- which(clu$csize == max(clu$csize)) # largest component(s)
  
  return(data.frame(NComp, CSize, nLarge = length(largclu)))
  
}