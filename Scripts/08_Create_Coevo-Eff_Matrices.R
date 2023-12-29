# --------- T-matrix (total effects) and Q-matrix (direct effects) --------

# Author: Kate P Maia
# Checked: 12/2023

# Computes and prints T and Q - the matrices of total (direct and indirect) and direct effects, using three values of the strength of interactions as a source of selective pressure (m). 

# --------------------- Loading library, code and data --------------------

library(tidyverse)
source("./Scripts/Functions/totaleffect.R")

dataset <- read.table("./Data/Project_Dataset.txt", header = T, sep = "\t", stringsAsFactors = F)

net_struct <- read.table("./Outputs/04_Net-Level_Struct.txt", header = T, stringsAsFactors = F)

dataset <- dataset[match(net_struct$ID, dataset$ID),]
all(dataset$ID == net_struct$ID)

net_struct$IntType <- dataset$IntType 

path <- "./Data/Matrices/"; net_names <- dir(path)

# REMOVING multLC & starLC
multLC <- c("Herb_CuevasReyes", "IZ_PA_DIM100", "IZ_PA_PAcamp", "M_AF_002_02", "M_AF_002_11")
starLC <- c("IZ_PA_DIM10", "IZ_PA_PA10", "M_AF_002_07", "M_AF_002_16", "M_PL_061_33")
net_names <- net_names[-grep(paste(c(multLC, starLC), collapse = "|"), net_names)]; length(net_names)

net_list <- lapply(paste(path, net_names, sep = ""), read.table, header = TRUE, sep = "\t", check.names = FALSE)

m_vec <- c(0.1, 0.5, 0.9)

#dir.create("./Outputs/TE-Matrices", showWarnings = TRUE)
#dir.create("./Outputs/Q-Matrices", showWarnings = TRUE)

# ---------------- Creates and prints binary TE ----------------

for (j in 1:length(m_vec)) {
  
  subf <- paste("Binary_m=", format(m_vec[j], nsmall = 1), sep = "")
  subf_TE <- paste("./Outputs/TE-Matrices/", subf, sep = "") # subfolder
  subf_Q <- paste("./Outputs/Q-Matrices/", subf, sep = "")
  #dir.create(subf_TE, showWarnings = TRUE)
  #dir.create(subf_Q, showWarnings = TRUE)
  
  for (i in 1:length(net_names)) {
    
    print(i)
    ID <- str_sub(net_names[i], end = -5)
    M <- net_list[[i]]; M[M > 0] <- 1 # all binary matrices
    
    TE <- totaleffect(M, m_vec[j])
    Q <- TE$Q; TE <- TE$TE 
    
    if (nrow(TE) != nrow(M) + ncol(M)) {print("ERROR IN TE DIMENSION")}
    if (any(rownames(TE) != c(rownames(M), colnames(M)))) {print("ERROR IN TE NAMES")}
    
    name <- paste(ID, "_TE_B_m=", m_vec[j], ".txt", sep = "")
    #write.table(TE, paste(subf_TE, "/", name, sep = ""), sep = "\t")
    name <- paste(ID, "_Q_B_m=", m_vec[j], ".txt", sep = "")
    #write.table(Q, paste(subf_Q, "/", name, sep = ""), sep = "\t")
  }
}