# -------------------- partition effect MdMt -------------------
# --------------------- nodepair_m6counter ---------------------

# Author: Kate P Maia
# Checked: 12/2023

# Partitions effects of node-pair-level dataframe (all within largest component after tij_dataframe.R and crosscomp_eff), into: inside module inside motif (inMdinMt), inside module between motifs (inMdbetMt), between modules inside motif (betMdinMt) and between modules between motifs (betMdbetMt). 
# input: a dataframe of node-pair-level effects (teff, deff, ieff), all within largest components, for all m.
# output: a summary dataframe with three rows (one per m) and five cols m, inMdinMt, inMdbetMt, betMdinMt, betMdbetMt which are the sum of effects in each partition.

#----------------------------------------------------------------

partition_eff_MdMt <- function(effects_df){
  
  m <- levels(factor(effects_df$m))
  result <- data.frame(m, inMdinMt = NA, inMdbetMt = NA, betMdinMt = NA, betMdbetMt = NA)
  
  for (j in 1:length(m)) {
    
    m_df <- effects_df[effects_df$m == m[j], c("OutSp", "InSp", "Effect", "InMod_Aff", "OutMod_Aff", "m6count")]
    
    inMd <- which(m_df$InMod_Aff == m_df$OutMod_Aff) # in and out in the same module
    betMd <- which(m_df$InMod_Aff != m_df$OutMod_Aff) # in and out in different modules
    inMt <- which(m_df$m6count != 0) # in and out in the same motif
    betMt <- which(m_df$m6count == 0) # in and out in different motifs
    
    test1 <- c(length(intersect(inMd, betMd)), length(intersect(inMt, betMt))) # any eff in & bet
    test2 <- c(length(inMd) + length(betMd), length(inMt) + length(betMt)) # all eff either in or bet
    if (any(test1 != 0)) {"ERROR: Partition 1"}
    if (any(test2 != nrow(m_df))) {"ERROR: Partition 2"}
    
    inMdinMt <- m_df[intersect(inMd, inMt),] # same module same motif
    inMdbetMt <- m_df[intersect(inMd, betMt),] # same module between motifs
    betMdinMt <- m_df[intersect(betMd, inMt),] # between modules same motif
    betMdbetMt <- m_df[intersect(betMd, betMt),] # between modules between motifs
    
    if (nrow(rbind(inMdinMt, inMdbetMt, betMdinMt, betMdbetMt)) != nrow(m_df)) {print("ERROR: Partition 3")}
    
    result[j, -1] <- c(sum(inMdinMt$Effect), sum(inMdbetMt$Effect), sum(betMdinMt$Effect), sum(betMdbetMt$Effect))
  }
  
  return(result)
}