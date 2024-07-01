# ---------- modsect_nm: summarises null model results ---------

# Author: Kate P Maia
# Checked: 12/2023

# Summarises null model results of module-sector congruence from dataframes into mean, sd and p-value of NM congruence values.
# input: a dataframe of empirical and NM module-sector data (congr_df).
# output: a summarised version of the input dataframe, in which columns with NM raw data are replaced with the mean, sd and p-value of null model results.

# --------------------------------------------------------------

modsect_nm <- function(df) {
  
  emp <- data.frame(congr = df[, 1]) # empirical 
  
  emp$mean <- apply(df[,-1], 1, mean) # mean and sd of NM values
  emp$sd <- apply(df[,-1], 1, sd)
  emp$pval <- colSums(apply(df, 1, function(x) x[-1] >= x[1])) / 1000
  
  ID <- rownames(emp); rownames(emp) <- NULL
  
  return(data.frame(ID, emp))
}
