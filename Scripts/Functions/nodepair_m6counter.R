# --------------------- nodepair_m6counter ---------------------

# Author: Kate P Maia
# Checked: 12/2023

# Counts the number of motif6 a node pair (rows in nodepair_df) participates in. For each nodepair_df row, it calls two auxiliary functions (see below).
# input: a dataframe with pairs of nodes extracted from teff/deff/ieff, and prf, the m6 profile of that network.
# output: nodepair_df* with the new column m6count.

# prf.long.f: Elongates the m6 profile of the network twice: first by changing it from a 4-column to a 2-column dataframe, second by adding a reversed version of each row (A-B, B-A) so that all pairs in nodepair_df are found regardless of species order. 

# aux1_m6counter: Matches each nodepair_df row with prf.long, to count how many times that pair (from nodepair_df) co-occurs in a motif. Called for each nodepair_df row.

#----------------------------------------------------------------

prf.long.f <- function(prf){
  colnames(prf) <- rep("name", 4)
  prf.long <- rbind(prf[, c(1,2)], prf[, c(1,3)], prf[, c(1,4)], prf[, c(2,3)], prf[, c(2,4)], prf[, c(3,4)])
  prf.long.rev <- prf.long[c(2,1)]
  colnames(prf.long) <- colnames(prf.long.rev) <- c("OutSp", "InSp") # to match nodepair_df
  prf.long <- rbind(prf.long, prf.long.rev)
  return(prf.long)
}

aux1_m6counter <- function(x, prf.long){
  nodepair_match <- suppressMessages(plyr::match_df(prf.long, data.frame(t(x))))
  m6count <- nrow(nodepair_match)
  return(m6count)
}

nodepair_m6counter <- function(nodepair_df, prf) {
  
  prf.long <- prf.long.f(prf) # turns prf (one motif per row) into pairwise df (all possible, inc rev)
  if (nrow(prf.long) != nrow(prf)*12) {print("ERROR IN PRF.LONG")}
  
  # counting pairs of two species
  m6count <- apply(nodepair_df, 1, aux1_m6counter, prf.long) # counts m6 each pair in nodepair_df is in
  
  # counting pairs of one species: indirect effects are also shared between a species and itself 
  posit_pairof1 <- which(nodepair_df$OutSp == nodepair_df$InSp) # position of 1-sp pairs in nodepair_df
  prf_string <- c(prf[,1], prf[,2], prf[,3], prf[,4])
  count_pairof1 <- table(prf_string)
  posit_pairof1_withcount <- match(names(count_pairof1), nodepair_df$OutSp[posit_pairof1])
  nodepair_df[posit_pairof1[posit_pairof1_withcount],] # from all pairof1, the ones with m6count
  #cbind(nodepair_df[posit_pairof1[posit_pairof1_withcount],], data.frame(count_pairof1))
  m6count[posit_pairof1[posit_pairof1_withcount]] <- count_pairof1
  
  # adding reverse pairs (BA) to nodepair_df to match AB and BA m6counts
  nodepair_df$m6count <- m6count
  if (any(nodepair_df$m6count[posit_pairof1[posit_pairof1_withcount]] != count_pairof1)) {print("ERROR in pair of 1")}
  nodepair_rev <- data.frame(OutSp = nodepair_df$InSp, InSp = nodepair_df$OutSp, m6count)
  nodepair_df <- rbind(nodepair_df, nodepair_rev)
  nodepair_df <- nodepair_df[!duplicated(nodepair_df),] # remove dup 1-sp pairs: AA is a single cell in an eff matrix
  rownames(nodepair_df) <- 1:nrow(nodepair_df)
  
  return(nodepair_df)
}