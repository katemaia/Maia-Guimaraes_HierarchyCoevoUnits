# ---------- modsect_nm: summarises null model results ---------

# Author: Kate P Maia
# Checked: 12/2023

# Summarises null model results of module-sector congruence from dataframes into mean, sd and p-value of NM congruence values.
# input: a dataframe of empirical and NM module-sector data (congr_df, congr_sp_df).
# output: a summarised version of the input dataframe, in which columns with NM raw data are replaced with the mean, sd and p-value of null model results.

# --------------------------------------------------------------

modsect_nm <- function(df) {
  
  emp <- df[, -which(colnames(df) %in% 1:1000)] # info + empirical 
  
  test <- df[, -c(1:(ncol(df) - 1001))] # test: emp + NM
  if (ncol(test) != 1001) {print("ERROR IN CROP")}
  
  emp$mean <- apply(test[,-1], 1, mean) # mean and sd of NM values
  emp$sd <- apply(test[,-1], 1, sd)
  emp$pval <- colSums(apply(test, 1, function(x) x[-1] >= x[1])) / 1000
  
  return(emp)
}
