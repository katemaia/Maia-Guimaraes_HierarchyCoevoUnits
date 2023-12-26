# ------------------------- fiedler_GC -------------------------

# Author: Kate P Maia
# Checked: 12/2023

# fiedler_GC: calculates the fiedler vector of networks'largest components. 
# output: list with one df for each largest component; with the name of each node in the component and its fiedler value

library(igraph)

# --------------------------------------------------------------

fiedler_GC <- function(M) {
  
  # Adjacency from incidence matrix
  nA <- nrow(M); nP <- ncol(M)
  q2 <- q4 <- M 
  q1 <- mat.or.vec(nA, nA); q3 <- mat.or.vec(nP, nP)
  upper <- cbind(q1, q2); lower <- cbind(t(q4), q3)
  colnames(upper) <- colnames(lower) <- c(rownames(M), colnames(M))
  A <- as.matrix(rbind(upper, lower))
  rownames(A) <- c(rownames(M), colnames(M))
  if (any(rownames(A) != colnames(A))) {print("ERROR IN FIEDLER_GC: 1")}
  
  g <- graph_from_adjacency_matrix(A, "undirected")
  clu <- igraph::components(g)
  largclu <- which(clu$csize == max(clu$csize)) # largest component(s)
  comp.sp <- groups(clu) # species affiliation to each component
  
  # For each largest component
  fiedler_GC <- list()
  for (i in 1:length(largclu)) {
    lc <- largclu[i] 
    lcsp <- comp.sp[[lc]]
    glc <- induced_subgraph(g, lcsp)
    if (any(V(glc)$name != lcsp)) {print("ERROR in glc graph")}
    laplaglc <- graph.laplacian(glc)
    eigen <- eigen(laplaglc)
    fied_position <- order(eigen$values)[2]
    fiedler <- eigen$vectors[, fied_position]
    fiedler_GC[[i]] <- data.frame(sp_GC = lcsp, fiedler = fiedler)
  }
  return(fiedler_GC)
}