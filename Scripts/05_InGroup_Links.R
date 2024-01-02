# --------- Net-level data: adds proportion of links inside groups --------

# Author: Kate P Maia
# Checked: 12/2023

# Calculates the proportion of links within (exchanged by species in the same) and between (exchanged by species in distinct) groups, for groups across network scales (sectors, modules and motifs). Adds results to net-level data. 

# --------------------- Loading library, code and data --------------------

library(tidyverse)
library(reshape2)
source("./Scripts/Functions/groupint_counter.R")
source("./Scripts/Functions/motifint_counter.R")

dataset <- read.table("./Data/Project_Dataset.txt", header = T, sep = "\t", stringsAsFactors = F)

net_struct <- read.table("./Outputs/04_Net-Level_Struct.txt", header = T, stringsAsFactors = F)
sp_struct <- read.table("./Outputs/04_Node-Level_Struct.txt", header = T, stringsAsFactors = F)

dataset <- dataset[match(net_struct$ID, dataset$ID),]
all(dataset$ID == net_struct$ID)

sp_struct <- sp_struct[order(sp_struct$ID),]
all(dataset$ID == unique(sp_struct$ID))

net_struct$IntType <- dataset$IntType 

# REMOVING multLC & starLC & TOO-LARGE (multiple largest components, star-shaped largest components)
multLC <- c("Herb_CuevasReyes", "IZ_PA_DIM100", "IZ_PA_PAcamp", "M_AF_002_02", "M_AF_002_11")
starLC <- c("IZ_PA_DIM10", "IZ_PA_PA10", "M_AF_002_07", "M_AF_002_16", "M_PL_061_33")
too_large <- "M_PL_062" # Removing very large network

net_names <- dir("./Data/Matrices")
net_names <- net_names[-grep(paste(c(multLC, starLC, too_large), collapse = "|"), net_names)]
matrices <- lapply(paste("./Data/Matrices/", net_names, sep = ""), read.table, header = TRUE, sep = "\t", check.names = "FALSE"); length(matrices)

prf_names <- dir("./Outputs/m6-Profiles")
prf_names <- grep(paste(net_names, collapse = "|"), prf_names, value = TRUE)
prf_list <- lapply(paste("./Outputs/m6-Profiles/", prf_names, sep = ""), read.table, header = TRUE, sep = "\t"); length(prf_list)

# ------ Computing interactions within and between groups ------

interaction_df <- data.frame(ID = character(0), Sct = numeric(0), Mod = numeric(0), Mot = numeric(0))

for (i in 1:length(net_names)){
  
  print(i)
  ID <- gsub(".txt", "",net_names[i])
  M <- matrices[[grep(ID, net_names)]]; prf <- prf_list[[grep(ID, prf_names)]]
  
  # largest component M - CHECKED
  sp_df <- sp_struct[sp_struct$ID == ID & !is.na(sp_struct$GC_Aff),]
  if (any(sp_df$GC_Aff != 1)) {print("ERROR")}
  rowsp <- sp_df[sp_df$MDim == 1, "Species"]; colsp <- sp_df[sp_df$MDim == 2, "Species"]
  M_GC <- M[rownames(M) %in% rowsp, colnames(M) %in% colsp]; M_GC[M_GC != 0] <- 1
  
  # Cross-Sector Interactions
  crosssct <- groupint_counter(sp_df, M_GC, "sector")
  Sct <- round(sum(crosssct[-1]) / crosssct[1], digits = 4)
  
  # Cross-Module Interactions
  crossmod <- groupint_counter(sp_df, M_GC, "module")
  Mod <- round(sum(crossmod[-1]) / crossmod[1], digits = 4)
  
  # Cross-Motif Interactions
  crossmot <- motifint_counter(prf, M_GC)
  Mot <- round(crossmot[-1] / crossmot[1], digits = 4)
  
  if (length(unique(c(crosssct[1], crossmod[1], crossmot[1]))) != 1) {"ERROR in total int"}
  
  interaction_df <- rbind(interaction_df, data.frame(ID, LinkSct = Sct, LinkMod = Mod, LinkMot = Mot))
}

head(interaction_df); dim(interaction_df)
net_struct <- left_join(net_struct, interaction_df, by = "ID")
net_struct <- net_struct %>% relocate(IntType, .after = ID)
head(net_struct); dim(net_struct)

#write.table(net_struct, "./Outputs/05_Net-Level_Struct.txt", sep = "\t")