# ----- Partition effects within largest components (module/subgraph) -----

# Author: Kate P Maia
# Checked: 12/2023

# Partitions three effect measures (total, direct and indirect effects)  - per network and per m value - in four groups: effects inside modules and inside motifs (inMdinMt), inside modules but between motifs (inMdbetMt), between modules but inside motifs (betMdinMt) and between modules and between motifs  (betMdbetMt). It needs to be run in parts as it takes a long time to run.
# Prints partition_df with the four partitions per network (ID), interaction type (IntType), effect measure (teff, deff or ieff), and m (0.1, 0.5 or 0.9). 

# --------------------- Loading library, code and data --------------------

library(plyr) # must be loaded before tidyverse
library(tidyverse)
library(reshape2)
source("./Scripts/Functions/tij_dataframe.R")
source("./Scripts/Functions/crosscomp_eff.R")
source("./Scripts/Functions/create_nodepair.R")
source("./Scripts/Functions/nodepair_m6counter.R")
source("./Scripts/Functions/partition_eff_MdMt.R")

dataset <- read.table("./Data/Project_Dataset.txt", header = T, sep = "\t", stringsAsFactors = F)

net_struct <- read.table("./Outputs/04_Net-Level_Struct.txt", header = T, stringsAsFactors = F)
sp_struct <- read.table("./Outputs/04_Node-Level_Struct.txt", header = T, stringsAsFactors = F)

dataset <- dataset[match(net_struct$ID, dataset$ID),]
all(dataset$ID == net_struct$ID)

sp_struct <- sp_struct[order(sp_struct$ID),]
all(dataset$ID == unique(sp_struct$ID))

net_struct$IntType <- dataset$IntType 

# -------------------------------------------------------------------------
# ------------- Partitioning effects in modules and subgraphs -------------

path <- "./Data/Matrices/"; net_names <- dir(path)

# REMOVING multLC & starLC & TOO-LARGE (multiple largest components, star-shaped largest components)
multLC <- c("Herb_CuevasReyes", "IZ_PA_DIM100", "IZ_PA_PAcamp", "M_AF_002_02", "M_AF_002_11")
starLC <- c("IZ_PA_DIM10", "IZ_PA_PA10", "M_AF_002_07", "M_AF_002_16", "M_PL_061_33")
too_large <- "M_PL_062" # Removing very large network
net_names <- net_names[-grep(paste(c(multLC, starLC, too_large), collapse = "|"), net_names)]

prf_names <- dir("./Outputs/m6-Profiles"); length(net_names) == length(prf_names)

# ---------- CROPPING TO RUN ----------

crop <- net_struct[!net_struct$ID %in% c(multLC, starLC, "M_PL_062"),]
dim(crop); any(is.na(crop$n_m6))

crop[crop$n_m6 >= 100000,] # 3 - MISSING Pearse&Alterm
#crop_names <- crop[crop$n_m6 >= 50000, "ID"][2]; length(crop_names) # ***** CROP *****

#crop[crop$n_m6 <= 100,] # 3 - MISSING Pearse&Alterm
crop_names <- crop[crop$n_m6 > 100000, "ID"]; length(crop_names) # ***** CROP *****

net_names <-  grep(paste(crop_names, collapse = "|"), net_names, value = TRUE) # ***** CROP *****
prf_names <-  grep(paste(crop_names, collapse = "|"), prf_names, value = TRUE) # ***** CROP *****

prf_list <- lapply(paste("./Outputs/m6-Profiles/", prf_names, sep = ""), read.table, header = TRUE, sep = "\t")

all(str_sub(net_names, end = -5) == str_sub(prf_names, start = 5, end = -5))

m_vec <- c(0.1, 0.5, 0.9)

pathTE <- "./Outputs/TE-Matrices/" # total effects
TEfolders <- paste(rep("Binary_m=", length(m_vec)), m_vec, "/", sep = "")
TEnames <- paste(rep("_TE_B_m=", length(m_vec)), m_vec, sep = "")

pathQ <- "./Outputs/Q-Matrices/" # direct effects
Qfolders <- paste(rep("Binary_m=", length(m_vec)), m_vec, "/", sep = "")
Qnames <- paste(rep("_Q_B_m=", length(m_vec)), m_vec, sep = "")

partition_df <- data.frame(ID = character(0), IntType = character(0), effmeasure = character(0), m = numeric(0), inMdinMt = numeric(0), inMdbetMt = numeric(0), betMdinMt = numeric(0), betMdbetMt = numeric(0))

# ***ATENTION: SLOW CODE. Add network subset here if needed

ptm <- proc.time()
for (i in 1:length(net_names)) {
  
  print(i)
  ID <- str_sub(net_names[i], end = -5)
  type <- net_struct[match(ID, net_struct$ID), "IntType"]
  prf <- prf_list[[i]]
  
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
  
  # Creates list of indirect effects: removes Q (ieffsum), then zeroes cells of direct effects (partition) 
  INDlist <- mapply(function(X, Y) {X[Y != 0] <- 0; return(X)}, X = TElist, Y = Qlist, SIMPLIFY = FALSE)
  if (any(c(nempty2 + nlink) != unlist(lapply(INDlist, function(x) sum(x == 0))))) {print("ERROR in IND")}
  
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
  
  # Create species pairs and counts how many times pairs co-occur in motifs 
  nodepair_df <- create_nodepair(teff, deff, ieff) # only AB pairs so counts for AB match counts for BA
  nodepair_df <- nodepair_m6counter(nodepair_df, prf)
  teff <- merge(teff, nodepair_df, by = c("OutSp", "InSp")) # merge: CHECKED
  deff <- merge(deff, nodepair_df, by = c("OutSp", "InSp")) 
  ieff <- merge(ieff, nodepair_df, by = c("OutSp", "InSp"))
  if (any(teff$m6count != deff$m6count | teff$m6count != ieff$m6count)) {print("ERROR IN M6Count")}
  
  # Partitions effects in inMdinMt, inMdbetMt, betMdinMt, betMdbetMt per m
  teff <- partition_eff_MdMt(teff); deff <- partition_eff_MdMt(deff); ieff <- partition_eff_MdMt(ieff)
  
  effects_df <- rbind(teff, deff, ieff)
  effects_df <- cbind(ID = ID, IntType = type, effmeasure = rep(c("teff", "deff", "ieff"), each = length(m_vec)), effects_df)
  
  partition_df <- rbind(partition_df, effects_df)
}
proc.time() - ptm

head(partition_df); dim(partition_df)
#write.table(partition_df, "./Outputs/10_Partition-Effects-MdMt-Dataframe.txt", sep = "\t")