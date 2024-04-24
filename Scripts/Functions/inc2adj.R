# ------------------------- inc2adj --------------------------

# Author: Kate P Maia
# Date: 03/2021

# Function which turns incidence into adjacency matrices, in which quadrants 1 and 3 are filled with zeros, and quadrants 2 and 4 correspond to the incidence matrix.
# input: M (incidence interaction matrix)
# output: A (adjacency interaction matrix).

# ------------------------------------------------------------

inc2adj <- function(M) {
  
  nA <- nrow(M); nP <- ncol(M)
  q1 <- mat.or.vec(nA, nA); q3 <- mat.or.vec(nP, nP)
  q2 <- q4 <- M 
  upper <- cbind(q1, q2); lower <- cbind(t(q4), q3)
  colnames(upper) <- colnames(lower) <- c(rownames(M), colnames(M))
  A <- as.matrix(rbind(upper, lower))
  rownames(A) <- c(rownames(M), colnames(M))
  return(A)
  
}
