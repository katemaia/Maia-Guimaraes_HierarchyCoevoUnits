# ------ Partition effects within largest components (sector/module) ------

# Author: Kate P Maia
# Checked: 12/2023

# Partitions three effect measures (total, direct and indirect effects)  - per network and per m value - in four groups: effects inside sectors and inside modules (inSinM), inside sectors but between modules (inSbetM), between sectors but inside modules (betSinM) and between sectors and between modules (betSbetM). Prints partition_df with the four partitions per network (ID), interaction type (IntType), effect measure (teff, deff or ieff), and m (0.1, 0.5 or 0.9).

# --------------------- Loading library, code and data --------------------

library(reshape2)
library(tidyverse)
source("./Scripts/Functions/tij_dataframe.R")
source("./Scripts/Functions/crosscomp_eff.R")
source("./Scripts/Functions/partition_eff.R")

dataset <- read.table("./Data/Project_Dataset.txt", header = T, sep = "\t", stringsAsFactors = F)

net_struct <- read.table("./Outputs/04_Net-Level_Struct.txt", header = T, stringsAsFactors = F)
sp_struct <- read.table("./Outputs/04_Node-Level_Struct.txt", header = T, stringsAsFactors = F)

dataset <- dataset[match(net_struct$ID, dataset$ID),]
all(dataset$ID == net_struct$ID)

sp_struct <- sp_struct[order(sp_struct$ID),]
all(dataset$ID == unique(sp_struct$ID))

net_struct$IntType <- dataset$IntType 

# -------------------------------------------------------------------------
# --- Partitioning effects in sectors and modules - DONE (SKIP TO L109) ---

path <- "./Data/Matrices/"; net_names <- dir(path)

# REMOVING multLC & starLC (see decide_N.R)
multLC <- c("Herb_CuevasReyes", "IZ_PA_DIM100", "IZ_PA_PAcamp", "M_AF_002_02", "M_AF_002_11")
starLC <- c("IZ_PA_DIM10", "IZ_PA_PA10", "M_AF_002_07", "M_AF_002_16", "M_PL_061_33")
net_names <- net_names[-grep(paste(c(multLC, starLC), collapse = "|"), net_names)]; length(net_names)

m_vec <- c(0.1, 0.5, 0.9)

pathTE <- "./Outputs/TE-Matrices/" # total effects
TEfolders <- paste(rep("Binary_m=", length(m_vec)), m_vec, "/", sep = "")
TEnames <- paste(rep("_TE_B_m=", length(m_vec)), m_vec, sep = "")

pathQ <- "./Outputs/Q-Matrices/" # direct effects
Qfolders <- paste(rep("Binary_m=", length(m_vec)), m_vec, "/", sep = "")
Qnames <- paste(rep("_Q_B_m=", length(m_vec)), m_vec, sep = "")

partition_df <- data.frame(ID = character(0), IntType = character(0), effmeasure = character(0), m = numeric(0), inSinM = numeric(0), inSbetM = numeric(0), betSinM = numeric(0), betSbetM = numeric(0))

for (i in 1:length(net_names)) {
  
  print(i)
  ID <- str_sub(net_names[i], end = -5)
  type <- net_struct[match(ID, net_struct$ID), "IntType"]
  
  # Reads total effects matrices - list of 3 (one per m value)
  TElist <- paste(rep(pathTE, 3), TEfolders, rep(ID, 3), TEnames, rep(".txt", 3), sep = "") %>% 
    lapply(read.table, header = TRUE, sep = "\t", check.names = FALSE) %>% lapply(as.matrix)
  
  # Reads direct effects matrices - list of 3 (one per m value)
  Qlist <- paste(rep(pathQ, 3), Qfolders, rep(ID, 3), Qnames, rep(".txt", 3), sep = "") %>% 
    lapply(read.table, header = TRUE, sep = "\t", check.names = FALSE) %>% lapply(as.matrix)
  
  # Confirms that Q effects match the number of links in M
  nS <- net_struct[net_struct$ID == ID, c("NRow", "NCol")]; nS <- nS$NRow + nS$NCol # nspecies
  nlink <- net_struct[net_struct$ID == ID, "NLinks"]*2 # nlinks in two quadrants 
  if (any(unlist(lapply(Qlist, function(x) {sum(x != 0)})) != nlink)) {print("ERROR in LINK")}
  
  # Removes diag from total effects matrix (removal that cascades to indirect effects matrix)
  nempty1 <- unlist(lapply(TElist, function(x) sum(x == 0))) # empty in TE
  TElist <- lapply(TElist, function(x) {diag(x) <- 0; return(x)})
  nempty2 <- unlist(lapply(TElist, function(x) sum(x == 0))) # empty in TE
  if (any(c(nempty1 + nS) != nempty2)) {print("ERROR in DIAG")}
  
  teffsum <- unlist(lapply(TElist, sum)); deffsum <- unlist(lapply(Qlist, sum)) # TE* and Q matrix sum
  
  # Creates list of indirect effects: removes Q (ieffsum), then zeroes cells of direct effects (partition) 
  INDlist <- mapply(function(X, Y) {X[Y != 0] <- 0; return(X)}, X = TElist, Y = Qlist, SIMPLIFY = FALSE)
  if (any(c(nempty2 + nlink) != unlist(lapply(INDlist, function(x) sum(x == 0))))) {print("ERROR in IND")}
  ieffsum <- unlist(lapply(INDlist, sum)) # IND* matrix sum
  
  # Starts to partition effects into groups: inin, inbet, betin, betbet
  sp_df <- sp_struct[sp_struct$ID == ID,] # CHECKED - matches dim of effect matrices
  m <- rep(m_vec, each = nrow(sp_df)^2)
  
  # Effects at the level of node pairs
  teff <- do.call(rbind, lapply(TElist, tij_dataframe, is.in = F, sp_df))
  deff <- do.call(rbind, lapply(Qlist, tij_dataframe, is.in = F, sp_df))
  ieff <- do.call(rbind, lapply(INDlist, tij_dataframe, is.in = F, sp_df))
  
  nrowdf <- nrow(TElist[[1]])^2*3
  if (any(c(nrow(teff), nrow(deff), nrow(ieff)) %% nrowdf != 0)) {print("ERROR-DF")}
  
  # Adds m column
  teff <- cbind(m, teff); deff <- cbind(m, deff); ieff <- cbind(m, ieff)
  
  # Tests whether cross-component effects exist and eliminates node-pairs not in largest component
  teff <- crosscomp_eff(teff); deff <- crosscomp_eff(deff); ieff <- crosscomp_eff(ieff)
  
  # Partitions effects in inSinM, inSbetM, betSinM, betSbetM per m
  teff <- partition_eff(teff); deff <- partition_eff(deff); ieff <- partition_eff(ieff)
  
  effects_df <- rbind(teff, deff, ieff); effects_df <- cbind(effects_df, sum = c(teffsum, deffsum, ieffsum))
  effects_df <- cbind(ID = ID, IntType = type, effmeasure = rep(c("teff", "deff", "ieff"), each = length(m_vec)), effects_df)
  
  partition_df <- rbind(partition_df, effects_df)
}

head(partition_df); dim(partition_df)
#write.table(partition_df, "./Outputs/09_Partition-Effects-SctMd-Dataframe.txt", sep = "\t")