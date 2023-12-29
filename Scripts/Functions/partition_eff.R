# ---------------------- partition effect ----------------------

# Author: Kate P Maia
# Checked: 12/2023

# Partitions effects of node-pair-level dataframe (all within largest component after tij_dataframe.R and crosscomp_eff), into: inside sector inside module (inSinM), inside sector between modules (inSbetM), between sectors inside module (betSinM) and between sectors between modules (betSbetM). 
# input: a dataframe of node-pair-level effects (teff, deff, ieff), all within largest components, for all m.
# output: a summary dataframe with three rows (one per m) and five cols m, inSinM, inSbetM, betSinM, betSbetM with the sum of  effects in each partition.

# --------------------------------------------------------------

partition_eff <- function(nodepair_df){
  
  m <- levels(factor(nodepair_df$m))
  result <- data.frame(m, inSinM = NA, inSbetM = NA, betSinM = NA, betSbetM = NA)
  
  for (j in 1:length(m)) {
    
    m_df <- nodepair_df[nodepair_df$m == m[j], -1]
    
    inS <- which(m_df$InFiedler == m_df$OutFiedler) # in and out in the same sector
    betS <- which(m_df$InFiedler != m_df$OutFiedler) # in and out in different sectors
    inM <- which(m_df$InMod_Aff == m_df$OutMod_Aff) # in and out in the same module
    betM <- which(m_df$InMod_Aff != m_df$OutMod_Aff) # in and out in different modules
    
    test1 <- c(length(intersect(inS, betS)), length(intersect(inM, betM))) # any eff in & bet
    test2 <- c(length(inS) + length(betS), length(inM) + length(betM)) # all eff either in or bet
    if (any(test1 != 0)) {"ERROR: Partition 1"}
    if (any(test2 != nrow(m_df))) {"ERROR: Partition 2"}
    
    inSinM <- m_df[intersect(inS, inM),] # same sector same module
    inSbetM <- m_df[intersect(inS, betM),] # same sector between modules
    betSinM <- m_df[intersect(betS, inM),] # between sectors same module
    betSbetM <- m_df[intersect(betS, betM),] # between sectors betwee modules
    
    if (nrow(rbind(inSinM, inSbetM, betSinM, betSbetM)) != nrow(m_df)) {print("ERROR: Partition 3")}
    
    result[j, -1] <- c(sum(inSinM$Effect), sum(inSbetM$Effect), sum(betSinM$Effect), sum(betSbetM$Effect))
  }
  
  return(result)
}
