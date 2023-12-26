# ------------ criteria_sp - connectivity parameter ------------

# Author: Kate P Maia
# Checked: 12/2023

# CGInt_count: computes.
# input: M (incidence interaction matrix [I, J])
# output: df with the name of each node in the network (sp), its fiedler value (Fiedler, see fiedler_GC.R), its GC affiliation if any (GC_aff: the number of the GC each node belongs to or NA if the node belongs to no GC), and the number of each node's GC cross-group interactions (CG_int: 0 if the node belongs to one GC but has no CG interactions and NA if the node belongs to no GC).

source("./Scripts/functions/fiedler_GC.R")

# --------------------------------------------------------------

# M <- read.table("./Data/Matrices/M_PL_001.txt") # matrix for tests

CGInt_count <- function(M) {
  
  M[M > 0] <- 1
  
  sp_df <- data.frame(sp = c(rownames(M), colnames(M)), Fiedler = NA, GC_aff = NA, CG_int = NA)
  
  int_df <- as.data.frame(which(M == 1, arr.ind = TRUE, useNames = FALSE)) # interaction dataframe
  int_df$V1 <- rownames(M)[int_df$V1]; int_df$V2 <- colnames(M)[int_df$V2]
  
  fiedler <- fiedler_GC(M) # list with GC dataframes (species and fiedler value)
  
  for (i in 1:length(fiedler)) { # For each GC in M
    
    # GC species dataframe
    sp_df_GC <- fiedler[[i]]
    sp_df_GC$fiedler <- ifelse(sp_df_GC$fiedler > 0, 1, -1)
    sp_df_GC$GC_aff <- i # species affiliation to GC
    
    sp_df[match(sp_df_GC$sp_GC, sp_df$sp), c("Fiedler", "GC_aff")] <- sp_df_GC[, c("fiedler", "GC_aff")]
  }
  return(sp_df)
}
