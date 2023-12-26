# --------------------------- motif6 ---------------------------

# Author: Kate P Maia
# Checked: 12/2023

# Lists all bipartite motif-6 (see Simmons et al. 2019) in M. Motif-6 is a fully-connected four species motif, being two species from each set.  
# output: profile, a df with four columns: the row and column number of the species in each motif

#----------------------------------------------------------------
 
motif6 <- function(M) {
  
  M[M != 0] <- 1
  motif1 <- as.data.frame(which(M != 0, arr.ind = TRUE)) # all motif1 (interacting species pair)
  rownames(motif1) <- NULL; colnames(motif1) <- c("R", "C")
  nmotif1 <- nrow(motif1)
  
  frame <- as.data.frame(t(combn(1:nmotif1, 2))) # on which all pairs of motif1 will be laid
  motif_df <- cbind(motif1[frame$V1,], motif1[frame$V2,]) # all possible comb of motif1
  motif_df <- as.data.frame(motif_df); colnames(motif_df) <- c("R1", "C1", "R2", "C2")
  
  # removes motifs with repeated sp: m6 requires 4sp
  requal <- which(motif_df$R1 == motif_df$R2)
  cequal <- which(motif_df$C1 == motif_df$C2)
  motif_df <- motif_df[-unique(c(requal, cequal)),]
  
  # creates profile by testing if 4sp in motif_df cross-interact (r1c2, r2c1)
  prfl <- data.frame(R1 = character(0), C1 = character(0), R2 = character(0), C2 = character(0))
  for (i in 1:nrow(motif_df)) {
    mi <- motif_df[i,]
    r1c2 <- motif1[motif1$R == mi$R1 & motif1$C == mi$C2,]
    r2c1 <- motif1[motif1$R == mi$R2 & motif1$C == mi$C1,]
    if (nrow(r1c2) == 1 & nrow(r2c1) == 1) {prfl <- rbind(prfl, mi)}
  } 
  prfl <- prfl[,c(1,3,2,4)]
  
  if (nrow(prfl) != 0) {
    prfl <- t(apply(prfl, 1, function(x) x <- cbind(x[1:2][order(x[1:2])], x[3:4][order(x[3:4])])))
    prfl <- prfl[!duplicated(prfl),]
    prfl <- as.data.frame(prfl)
    if (ncol(prfl) == 1) {prfl <- as.data.frame(t(prfl))}
    rownames(prfl) <- NULL; colnames(prfl) <- c("R1", "R2", "C1", "C2")
    return(prfl)
  } else {return(prfl)}
}