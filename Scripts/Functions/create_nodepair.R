# ----------------------- create_nodepair ----------------------

# Author: Kate P Maia
# Checked: 12/2023

# Performs three functions: i) checks if species pairs in teff, deff and ieff match, ii) checks if species pairs within teff, deff and ieff match across m values, and iii) returns a dataframe with unique node pairs. 
# input: teff, deff and ieff.
# output: a summary dataframe one node pair per row (only names).

# --------------------------------------------------------------

create_nodepair <- function(teff, deff, ieff) {
  
  td <- any(teff[, c("m", "OutSp", "InSp")] != deff[, c("m", "OutSp", "InSp")])
  ti <- any(teff[, c("m", "OutSp", "InSp")] != ieff[, c("m", "OutSp", "InSp")])
  
  if (td | ti) {print("ERROR ACROSS EFF DATAFRAMES")}
  
  tm <- nrow(unique(teff[,c("OutSp", "InSp")])) != nrow(teff)/3
  dm <- nrow(unique(deff[,c("OutSp", "InSp")])) != nrow(deff)/3
  im <- nrow(unique(ieff[,c("OutSp", "InSp")])) != nrow(ieff)/3
  
  if (tm | dm | im) {print("ERROR INSIDE EFF DATAFRAMES")}
  
  nodepair_df <- unique(teff[,c("OutSp", "InSp")])
  nodepair_df$OutSp <- as.character(nodepair_df$OutSp)
  nodepair_df$InSp <- as.character(nodepair_df$InSp)
  
  # removes reverse duplicates: order does not matter, A-B count matches B-A count
  nodepair_df <- nodepair_df[!duplicated(apply(nodepair_df, 1, function(x) paste(sort(x), collapse = ''))),] # sort in alphabetic order, collapses and removes duplicates
  rownames(nodepair_df) <- 1:nrow(nodepair_df)
  return(nodepair_df)                        
                        
}