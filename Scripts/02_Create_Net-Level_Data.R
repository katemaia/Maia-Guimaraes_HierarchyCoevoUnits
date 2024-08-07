# ------------------------ Network-Level Structure ------------------------

# Author: Kate P Maia
# Checked: 08/2024

# Creates net_struct dataframe containing information on the structure of cohesive species groups across ecological network scales. Divided in:
# Loading library, code and data
# COMPONENTS IN NETWORKS
# SECTORS IN LARGEST COMPONENTS
# MODULES IN LARGEST COMPONENTS
# MODULE-SECTOR CONGRUENCE

# output: df with network ID, number of rows, columns and links (NRow, NCol, NLinks), connectance, number of components (NComp), size of the largest component in number of species (LargeCompSize), number of largest components, number of modules (_NMod), modularity (_Modularity) and null model results of modularity analysis ran with Simulated Annealing (SA) and Fast Grid (FG) algorithms, and module-sector (MS) observed congruence (congr), and mean, standard deviation (sd) and p-value (pval) of null model.

# --------------------- Loading library, code and data --------------------

library(tidyverse)
library(GoodmanKruskal)
source("./Scripts/functions/component_net.R")
source("./Scripts/functions/modsect_nm.R")

dataset <- read.table("./Data/Project_Dataset.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
dataset$IntType <- factor(dataset$IntType) # network information

sp_struct <- read.table("./Outputs/01_Node-Level_Struct.txt", header = TRUE, sep = "\t")

# modularity results calculated with MODULAR
SA <- read.table("./Modularity/Matrices/SimulatedAnnealing/OUT_MOD.txt", header = TRUE) # no null model
FG <- read.table("./Modularity/Matrices/FastGrid/OUT_MOD.txt", header = TRUE)
FGNull <- read.table("./Modularity/Null_Matrices/FastGrid/OUT_MOD.txt", header = TRUE) # null networks

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

FGNull$File <- str_sub(FGNull$File, start = 3, end = -5) 
FGNull$File <- unlist(lapply(str_split(FGNull$File, "_GC1"), function(x) x[1]))
all(unique(FGNull$File) %in% SA$File); all(SA$File %in% unique(FGNull$File))

# Calculating FG null model results
FGNull <- lapply(split(FGNull, f = FGNull$File), function(x) x$Modularity)
FGEmp <- lapply(split(FG, f = FG$File), function(x) x$Modularity)
all(names(FGNull) == names(FGEmp))
FG_PNull <- mapply(function(X, Y) {sum(Y >= X)}, X = FGEmp, Y = FGNull) / 100

# Adding modularity results to net_struct
colnames(SA)[1:3] <- colnames(FG)[1:3] <- c("ID", "NMod", "Modularity")
colnames(SA)[2:3] <- paste0("SA_", colnames(SA)[2:3])
colnames(FG)[2:3] <- paste0("FG_", colnames(FG)[2:3])

net_struct <- left_join(net_struct, SA[,1:3], by = "ID")
net_struct <- left_join(net_struct, FG[,1:3], by = "ID")
net_struct$FG_PNull <- NA
net_struct$FG_PNull <- FG_PNull[match(net_struct$ID, names(FG_PNull))]

head(net_struct); dim(net_struct)

# ------------------------ MODULE-SECTOR CONGRUENCE -----------------------

sp_struct <- sp_struct[!sp_struct$ID %in% c(multLC, starLC),]; length(unique(sp_struct$ID))
sp_struct <- sp_struct[!is.na(sp_struct$Mod_Aff),] # only sp in LC

sp_list <- split(sp_struct, f = sp_struct$ID) # Empirical mod-sect congruence
sp_list <- lapply(sp_list, function(x) x[, c("Fiedler", "Mod_Aff")]) # sp aff to sector and module
congr <- do.call(rbind, lapply(sp_list, function(x) GKtau(x$Mod_Aff, x$Fiedler)))

NM <- matrix(NA, nrow(congr), 1000) # Adding null model; stores emp value in col 1
congr_df <- cbind(congr$tauxy, NM); rownames(congr_df) <- rownames(congr)

for (i in 1:nrow(congr_df)) { # aprox 5'
  
  print(i)
  sp_M <- sp_struct[sp_struct$ID == rownames(congr_df)[i],] # component sp 
  if (any(is.na(sp_M)) | any(sp_M$GC_Aff != 1)) {print("ERROR IN sp_M")}
  
  for (j in 1:1000) { # null model
    
    sp_M$Fiedler <- unlist(tapply(sp_M$Fiedler, sp_M$MDim, sample))
    tau <- GKtau(sp_M$Mod_Aff, sp_M$Fiedler)
    congr_df[i, j + 1] <- tau$tauxy
    
  }
} 

congr_df <- modsect_nm(congr_df); colnames(congr_df)[-1] <- paste0("MS_", colnames(congr_df)[-1])

net_struct <- left_join(net_struct, congr_df, by = "ID")
head(net_struct); dim(net_struct)
#write.table(net_struct, "./Outputs/02_Net-Level_Struct.txt", sep = "\t")