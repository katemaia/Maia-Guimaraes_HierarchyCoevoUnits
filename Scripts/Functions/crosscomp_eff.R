# ------------------ cross component effects -------------------

# Author: Kate P Maia
# Checked: 12/2023

# Tests whether a node-pair-level dataframe (output of tij_dataframe.R) contains effects between components (cross component effects do not exist). If so, prints error message, otherwise removes cross-component node pairs and prints the dataframe.
# input: a dataframe of node-pair-level effects (tout, tin, dout, din, iout, iin).
# output: the input dataframe, but without rows (node-pairs) in which nodes belong to different components (all from GC_Aff = 1), and without node's component affiliation columns (GC_Aff),

# --------------------------------------------------------------

crosscomp_eff <- function(nodepair_df){
  
  inout <- any(nodepair_df[is.na(nodepair_df$InGC_Aff) & !is.na(nodepair_df$OutGC_Aff), "Effect"] != 0)
  outin <- any(nodepair_df[is.na(nodepair_df$OutGC_Aff) & !is.na(nodepair_df$InGC_Aff), "Effect"] != 0)
  if (any(inout, outin)) {print("ERROR: Effects between components")}
  
  nodepair_df <- nodepair_df[!is.na(nodepair_df$InGC_Aff) & !is.na(nodepair_df$OutGC_Aff),] # out of GC
  nodepair_df <- nodepair_df[nodepair_df$InGC_Aff == 1 & nodepair_df$OutGC_Aff == 1,]
  if (any(c(nodepair_df$InGC_Aff, nodepair_df$OutGC_Aff) != 1)) {print("ERROR: Other components")}
  
  return(nodepair_df[, -c(6, 10)])
  
}