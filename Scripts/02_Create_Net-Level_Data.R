# ------------------------ Network-Level Structure ------------------------

# Author: Kate P Maia
# Checked: 12/2023

# Creates net_struct dataframe containing information on the structure of cohesive species groups across ecological network scales. Divided in:
# Loading library, code and data
# COMPONENTS IN NETWORKS
# SECTORS IN LARGEST COMPONENTS
# MODULES IN LARGEST COMPONENTS
# MODULE-SECTOR CONGRUENCE
# SUBGRAPHS IN LARGEST COMPONENTS

# output: df with network ID, number of rows, columns and links (NRow, NCol, NLinks), connectance, number of components (NComp), size of the largest component in number of species (LargeCompSize), number of largest components, number of modules (_NMod), modularity (_Modularity) and null model results of modularity analysis ran with Simulated Annealing (SA) and Fast Grid (FG) algorithms, and module-sector (MS) observed congruence (congr), and mean, standard deviation (sd) and p-value (pval) of null model.

# --------------------- Loading library, code and data --------------------

library(tidyverse)
source("./Scripts/functions/component_net.R")
source("./Scripts/functions/modsect_nm.R")

dataset <- read.table("./Data/Project_Dataset.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
dataset$IntType <- factor(dataset$IntType) # network information

sp_struct <- read.table("./Outputs/01_Node-Level_Struct.txt", header = TRUE, sep = "\t")

# modularity results calculated with MODULAR
SA <- read.table("./Modularity/Matrices/SimulatedAnnealing/OUT_MOD.txt", header = TRUE) # no null model
FG <- read.table("./Modularity/Matrices/FastGrid/OUT_MOD.txt", header = TRUE) # with null model

path <- "./Data/Matrices/"; net_names <- dir(path) # matrices
net_list <- lapply(paste(path, net_names, sep = ""), read.table, header = TRUE, sep = "\t", check.names = FALSE)

# ------------------------- COMPONENTS IN NETWORKS ------------------------

net_struct <- data.frame(ID = character(0), NRow = numeric(0), NCol = numeric(0), NLinks = numeric(0), Connectance = numeric(0), NComp = numeric(0), LargeCompSize = numeric(0), nLargeComp = numeric(0))

for (i in 1:length(net_list)) {
  
  M <- net_list[[i]]; M[M != 0] <- 1
  ID <- str_sub(net_names[i], end = -5)
  NRow <- nrow(M)
  NCol <- ncol(M)
  NLinks <- sum(M)
  Conn <- round(NLinks / (NRow * NCol), digits = 4)
  
  cmpnt <- component_net(M)
  df <- data.frame(ID = ID, NRow = NRow, NCol = NCol, NLinks = NLinks, Connectance = Conn, NComp = cmpnt$NComp, LargeCompSize = cmpnt$CSize, nLargeComp = cmpnt$nLarge)
  net_struct <- rbind(net_struct, df)
  
}

head(net_struct); dim(net_struct)

# --------------------- SECTORS IN LARGEST COMPONENTS ---------------------

# No measurement at network-level, see 02_Create_Node-Level_Data.R.

# --------------------- MODULES IN LARGEST COMPONENTS ---------------------

# REMOVING multLC (with multiple largest components) & starLC (with star-shaped largest components)
multLC <- c("Herb_CuevasReyes", "IZ_PA_DIM100", "IZ_PA_PAcamp", "M_AF_002_02", "M_AF_002_11")
starLC <- c("IZ_PA_DIM10", "IZ_PA_PA10", "M_AF_002_07", "M_AF_002_16", "M_PL_061_33")
ID_order <- net_struct[!net_struct$ID %in% c(multLC, starLC), "ID"]

SA$File <- str_sub(SA$File, start = 3, end = -9) # cropping net ID to match net_struct
FG$File <- str_sub(FG$File, start = 3, end = -9) 

SA <- SA[match(ID_order, SA$File),]; FG <- FG[match(ID_order, FG$File),]

# Adding modularity results to net_struct
colnames(SA)[-4] <- colnames(FG)[-4] <- c("ID", "NMod", "Modularity", "PNull")
colnames(SA)[-1] <- paste0("SA_", colnames(SA)[-1])
colnames(FG)[-1] <- paste0("FG_", colnames(FG)[-1])

net_struct <- left_join(net_struct, SA[,1:3], by = "ID")
net_struct <- left_join(net_struct, FG[,c(1:3,5)], by = "ID")

head(net_struct); dim(net_struct)

# ------------------------ MODULE-SECTOR CONGRUENCE -----------------------

sp_struct <- sp_struct[!sp_struct$ID %in% c(multLC, starLC),]; length(unique(sp_struct$ID))
sp_struct <- sp_struct[!is.na(sp_struct$Mod_Aff),] # only sp in LC

# Btw 0.5 (all modules perfectly split) and 1 (all modules fully contained in dominant sector).
congr_df <- sp_struct
congr_df <- unique(congr_df[, c("ID", "GC_Aff")]); dim(unique(congr_df))
congr_df$congr <- NA

NM <- matrix(NA, nrow(congr_df), 1000); congr_df <- cbind(congr_df, NM) # Adding null model

for (i in 1:nrow(congr_df)) { # around 10'
  
  print(i)
  sp_M <- sp_struct[sp_struct$ID == congr_df$ID[i] & sp_struct$GC_Aff == congr_df$GC_Aff[i],] # comppnent sp 
  if (any(is.na(sp_M))) {print("ERROR IN sp_M")}
  df <- as.data.frame.matrix(table(sp_M$Mod_Aff, sp_M$Fiedler)) # N sp in [module, sector]
  df <- apply(df/rowSums(df), 1, max) # max congruence with sector (in %) per module
  congr_df[i, "congr"] <- mean(df)

  for (j in 1:1000) { # null model
    
    sp_M$Fiedler <- unlist(tapply(sp_M$Fiedler, sp_M$MDim, sample))
    df <- as.data.frame.matrix(table(sp_M$Mod_Aff, sp_M$Fiedler)) # N sp in [module, sector]
    df <- apply(df/rowSums(df), 1, max) # max congruence with sector (in %) per module
    congr_df[i, which(colnames(congr_df) == j)] <- mean(df)
    
  }
} 

congr_df <- modsect_nm(congr_df); head(congr_df) # summary NM results
colnames(congr_df)[3:6] <- paste0("MS_", colnames(congr_df)[3:6])

net_struct <- left_join(net_struct, congr_df[,-2], by = "ID")
head(net_struct); dim(net_struct)
#write.table(net_struct, "./Outputs/02_Net-Level_Struct.txt", sep = "\t", row.names = FALSE)