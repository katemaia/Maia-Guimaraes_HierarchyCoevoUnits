# ---------------------- Null networks for modularity ---------------------

# Author: Kate P Maia
# Checked: 06/2024

# Creates null networks to compute the level of modularity of empirical networks. Divided in:
# Loading library, code and data
# Running null models and saving null networks

# --------------------- Loading library, code and data --------------------

library(tidyverse)
library(bipartite)

dataset <- read.table("./Data/Project_Dataset.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# REMOVING multLC (with multiple largest components) & starLC (with star-shaped largest components)
multLC <- c("Herb_CuevasReyes", "IZ_PA_DIM100", "IZ_PA_PAcamp", "M_AF_002_02", "M_AF_002_11")
starLC <- c("IZ_PA_DIM10", "IZ_PA_PA10", "M_AF_002_07", "M_AF_002_16", "M_PL_061_33")

netID <- dataset$ID[!dataset$ID %in% c(multLC, starLC)]

path <- "./Modularity/Matrices/"
m.names <- grep(paste(netID, collapse = "|"), dir(path), value = TRUE)
m.list <- lapply(paste0(path, m.names), read.table, header = FALSE, sep = "\t")

# -------------- Running null models and saving null networks -------------

path <- "./Modularity/Null_Matrices/"

for (i in 1:length(m.list)) {
  
  print(i)
  ID <- gsub(".txt", "", m.names[i])
  memp <- m.list[[i]]
  
  null <- vaznull(100, memp)
  allconn <- unlist(lapply(null, function(x) all(rowSums(x) != 0) & all(colSums(x) != 0)))
  print(sum(allconn))
  
  #mapply(FUN = write.table, null, paste0(path, ID, "_", 1:length(null), ".txt"), row.names = FALSE, col.names = FALSE, sep = "\t", SIMPLIFY = FALSE)
  
}
