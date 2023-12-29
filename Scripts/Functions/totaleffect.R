# ------------------------ totaleffect -------------------------

# Author: Kate P Maia
# Checked: 12/2023

# Creates the matrix of total direct and indirect effects, following Pires et al. 2020. 
# input: M (interaction incidence matrix) and m (coefficient of decay).
# output: a list of matrices Q (adjacency direct effects matrix) and TE (adjacency total effects matrix).

# --------------------------------------------------------------

totaleffect <- function(M, m) {
  
  # Turns M into A (incidence into adjacency)
  Q2 <- M; Q4 <- t(M)
  if (is.null(colnames(Q4))) {colnames(Q4) <- rownames(Q2)}
  Q1 <- matrix(0, nrow(M), nrow(M)); rownames(Q1) <- colnames(Q1) <- rownames(M)
  Q3 <- matrix(0, ncol(M), ncol(M)); rownames(Q3) <- colnames(Q3) <- colnames(M)
  A <- rbind(cbind(Q1, Q2), cbind(Q4, Q3))
  if (any(rownames(A) != c(rownames(M), colnames(M)))) {"ERROR 1 in A"}
  if (any(colnames(A) != c(rownames(M), colnames(M)))) {"ERROR 2 in A"}
  
  D <- as.matrix(A/rowSums(A)); if (any(rowSums(D) != 1)) {"ERROR in D"} # D (dependence)
  
  R <- diag(rep(m, nrow(A))) # Creates R (coefficient of decay matrix)
  
  Q <- R %*% D; if (any(!near(rowSums(Q), m))) {"ERROR in Q"} # Computes Q (rescaled dependence)
  
  I <- diag(1, nrow(A)) # Creates I (identity)
  
  TE <- solve(I - Q) # Creates TE (total effect)
  colnames(TE) <- rownames(TE); rownames(Q) <- colnames(Q)
  if (any(rownames(TE) != rownames(A))) {print("ERROR in TE")}
  if (any(rownames(Q) != rownames(A))) {print("ERROR in Q")}
  
  return(list(Q = Q, TE = TE))
}