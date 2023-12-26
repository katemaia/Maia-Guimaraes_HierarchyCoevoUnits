# - component_sp - species affiliation to component and sector -

# Author: Kate P Maia
# Checked: 12/2023

# Computes species affiliation to networks largest components and sectors (using fiedler_vec.R).
# input: M (incidence interaction matrix [I, J])
# output: df with the name of each node in the network (sp), its largest affiliation if any (GC_aff: the number of the GC each node belongs to or NA if the node belongs to no GC), its fiedler value (Fiedler).

source("./Scripts/functions/fiedler_vec.R")

# --------------------------------------------------------------

component_sp <- function(M) {
  
  M[M > 0] <- 1
  
  sp_df <- data.frame(sp = c(rownames(M), colnames(M)), GC_aff = NA, Fiedler = NA)
  
  int_df <- as.data.frame(which(M == 1, arr.ind = TRUE, useNames = FALSE)) # interaction dataframe
  int_df$V1 <- rownames(M)[int_df$V1]; int_df$V2 <- colnames(M)[int_df$V2]
  
  fiedler <- fiedler_vec(M) # list with largest component dfs (species and fiedler value)
  
  for (i in 1:length(fiedler)) { # for networks with multiple largest components
    
    sp_df_GC <- fiedler[[i]] # GC species dataframe
    sp_df_GC$GC_aff <- i # species affiliation to largest component
    sp_df_GC$fiedler <- ifelse(sp_df_GC$fiedler > 0, 1, -1) # species affiliation to sector
    
    sp_df[match(sp_df_GC$sp_GC, sp_df$sp), c("GC_aff", "Fiedler")] <- sp_df_GC[, c("GC_aff", "fiedler")]
  }
  return(sp_df)
}
