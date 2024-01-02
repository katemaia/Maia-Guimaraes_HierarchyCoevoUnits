# ---------------------- groupint_counter ----------------------

# Author: Kate P Maia
# Checked: 12/2023

# Counts the number of links inside the network's largest component, and inside groups at other levels (sectors or modules). 
# input: sp_df (species level data for species in the largest component of a network), M_GC (interaction matrix of the network's largest component), group (sectors or modules).
# output: a vector with the number of interactions inside the largest component (All) and inside each group (G).

# --------------------------------------------------------------

groupint_counter <- function(sp_df, M_GC, group = c("sector", "module")){
  
  df <- sp_df[, c("Species", "MDim", ifelse(group == "sector", "Fiedler", "Mod_Aff"))]
  colnames(df)[3] <- "Group"
  groupint <- sum(M_GC) # all interactions
  
  # for each group at the selected level
  for (j in 1:length(unique(df$Group))){
    
    rowsp <- df[df$Group == unique(df$Group)[j] & df$MDim == 1,] # row species
    colsp <- df[df$Group == unique(df$Group)[j] & df$MDim == 2,] # col species 
    
    mgroup <- M_GC # crops M of the GC to become the M of the group
    mgroup <- mgroup[which(rownames(mgroup) %in% rowsp$Species),]
    mgroup <- mgroup[, which(colnames(mgroup) %in% colsp$Species)]
    
    dimmg <- dim(as.matrix(mgroup))
    if (dimmg[1] != nrow(rowsp) | dimmg[2] != nrow(colsp)) {print("ERROR GROUP DIM")}
    groupint <- c(groupint, sum(as.matrix(mgroup)))
    
  }
  
  names(groupint) <- c("All", paste("G", 0:(length(groupint) - 2), sep = ""))
  return(groupint)
}