# ---------------------------- Subgraph profile ---------------------------

# Author: Kate P Maia
# Checked: 12/2023

# Creates and prints the profile of motif-6 within network's largest components. Motif-6 (Simmons et al. 2019) corresponds to the subgraph selected to represent the smallest cohesive group in ecological networks.

# -------------------- Loading library, code and data --------------------

library(stringr)
source("./Scripts/Functions/motif6.R")

dataset <- read.table("./Data/Project_Dataset.txt", header = T, sep = "\t", stringsAsFactors = F)

sp_struct <- read.table("./Outputs/01_Node-Level_Struct.txt", header = T, stringsAsFactors = F)
net_struct <- read.table("./Outputs/02_Net-Level_Struct.txt", header = T, stringsAsFactors = F)

dataset <- dataset[match(net_struct$ID, dataset$ID),]
all(dataset$ID == net_struct$ID)

sp_struct <- sp_struct[order(sp_struct$ID),]
all(dataset$ID == levels(sp_struct$ID))

net_struct$IntType <- dataset$IntType 
sp_struct$IntType <- dataset[match(sp_struct$ID, dataset$ID), "IntType"]

path <- "./Data/Matrices/"; net_names <- dir(path)

# REMOVING multLC & starLC & TOO-LARGE (multiple largest components, star-shaped largest components)
multLC <- c("Herb_CuevasReyes", "IZ_PA_DIM100", "IZ_PA_PAcamp", "M_AF_002_02", "M_AF_002_11")
starLC <- c("IZ_PA_DIM10", "IZ_PA_PA10", "M_AF_002_07", "M_AF_002_16", "M_PL_061_33")
too_large <- "M_PL_062" # it was not possible to profile this matrix
net_names <- net_names[-grep(paste(c(multLC, starLC, too_large), collapse = "|"), net_names)] # N = 365

net_list <- lapply(paste(path, net_names, sep = ""), read.table, header = TRUE, sep = "\t", check.names = FALSE)

#dir.create("./Outputs/m6-Profiles", showWarnings = TRUE)

#------------- Creates and prints m6 profiles in networks' GCs --------------

for (i in 1:length(net_names)) {
  
  print(i)
  ID <- str_sub(net_names[i], end = -5)
  M <- net_list[[i]]
  
  sp_M <- sp_struct[sp_struct$ID == ID, c("ID", "Species", "MDim", "Fiedler", "GC_Aff", "Mod_Aff")]
  sp_M <- sp_M[!is.na(sp_M$GC_Aff) & sp_M$GC_Aff == 1,]
  
  M <- M[rownames(M) %in% sp_M$Species, colnames(M) %in% sp_M$Species]
  if (any(c(rownames(M), colnames(M)) != sp_M$Species)) {print("ERROR")}
  
  profile <- motif6(M) # all motifs in GC M
  profile_names <- profile
  
  profile_names$R1 <- rownames(M)[profile$R1]
  profile_names$R2 <- rownames(M)[profile$R2]
  profile_names$C1 <- colnames(M)[profile$C1]
  profile_names$C2 <- colnames(M)[profile$C2]
  
  #write.table(profile_names, paste("./Outputs/m6-Profiles/prf_", ID, ".txt", sep = ""), sep = "\t")
}
