# ------------------------ Tij dataframe -----------------------

# Author: Kate P Maia
# Checked: 12/2023

# Creates dataframe with the effects between nodes i and j. Effects are total (TE), direct (Q) or indirect (IND) depending on the input effects matrix, and in or out depending on is.in argument (in: is.in = T or out: is.in = F).  
# input: an effects matrix, and is.in (T or F).
# output: a dataframe where each row is a species pair with the effect between them, and the dimension (MDim), component (GC_Aff), sector (Fiedler) and module (Mod_Aff) affiliation of each species in the pair. 
# If is.in = F, matrix is transposed twice, and out effects are calculated (col 1 outsp, col 2 insp).
# If is.in = T, matrix is transposed once, and in effects are calculated (col 1 insp, col 2 outsp).

# --------------------------------------------------------------

tij_dataframe <- function(matrix, is.in = FALSE, sp_df){
  
  # trick so melt uses matrix rownames 
  if (is.in) {matrix <- t(matrix)} else {matrix <- t(t(matrix))}
  
  df <- melt(matrix); df <- df[, c(2, 1, 3)]
  if (is.in) {colnames(df) <- c("InSp", "OutSp", "Effect")} else {
    colnames(df) <- c("OutSp", "InSp", "Effect")}
  
  insp_aff <- sp_df[match(df$InSp, sp_df$Species), c("MDim", "GC_Aff", "Fiedler", "Mod_Aff")]
  outsp_aff <- sp_df[match(df$OutSp, sp_df$Species), c("MDim", "GC_Aff", "Fiedler", "Mod_Aff")]
  colnames(insp_aff) <- paste("In", colnames(insp_aff), sep = "")
  colnames(outsp_aff) <- paste("Out", colnames(outsp_aff), sep = "")
  
  df <- cbind(df, insp_aff, outsp_aff)
  return(df)
}