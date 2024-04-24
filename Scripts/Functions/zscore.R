# -------------------------- zscore --------------------------

# Author: Kate P Maia
# Date: 03/2024

# Function which computes the zscore and normalised zscore (Milo 2004) of congruence measures relative to a null model. zscores allow the comparison of observed congruence values with comparison to a null model, and normalisation facilitates the comparison of zscored values across networks of different sizes. 
# input: a dataframe with columns: network (ID), interaction type (Code), congruence values of empirical networks (Obs) and null model results (mean and sd of congruence values produced by the null model).
# output: a dataframe with columns: network (ID), interaction type (Code), zscore (Z) and normalised zscore (NormZ) values.

# ------------------------------------------------------------

zscore <- function(data) {
  
  colnames(data) <- c("ID", "Code", "Obs", "Mean", "Sd")
  data <- data %>% 
    mutate(Z = (Obs - Mean) / Sd) %>% # zscore
    filter(!is.infinite(Z)) %>% # Inf from 1/0
    mutate(Z2 = Z ^ 2) %>% # creates denominator for normalisation: squares, sums, sqrt 
    mutate(SumZ2 = sum(Z2, na.rm = TRUE)) %>% 
    mutate(SqrtSumZ2 = sqrt(SumZ2)) %>% 
    mutate(NormZ = Z / SqrtSumZ2) %>% 
    select(ID, Code, Z, NormZ)
  
  return(data)
  
}
