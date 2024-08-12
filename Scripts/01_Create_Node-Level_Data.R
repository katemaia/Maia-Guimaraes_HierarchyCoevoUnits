# -------------------------- Node-Level Structure -------------------------

# Author: Kate P Maia
# Checked: 08/2024

# Creates sp_struct dataframe containing information on the role of each species to the cohesive groups across ecological network scales, and prints matrices of networks' largest components for modularity analysis. Divided in:
# Loading library, code and data
# NODE AFFILIATION TO COMPONENTS AND SECTORS (prints LC matrices)
# NODE AFFILIATION TO MODULES

# output: df with network network ID, species' name, species' MDim (matrix dimension to identify whether a consumer or resource species), species' degree, largest component affiliation (GC_Aff), sector affiliation (Fiedler), and module affiliation (Mod_Aff).

# We replaced an igraph function which was deprecated in the auxiliary function fiedler_vec. The new function equally splits species in sectors, but inverts the sign of sector tags for some networks. We emphasize that the sign is completely arbitrary, therefore, not affecting paper results. To produce sector affiliation tags with equal signs as the ones used in the paper, uncomment lines 35 and 48 of this script.   

# --------------------- Loading library, code and data --------------------

library(tidyverse)
library(igraph)

source("./Scripts/functions/component_sp.R")

dataset <- read.table("./Data/Project_Dataset.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE) # dataset to provide interaction type

path <- "./Data/Matrices/"; net_names <- dir(path) # interaction matrices
net_names <- net_names[match(dataset$ID, str_sub(net_names, end = -5))]
all(dataset$ID == str_sub(net_names, end = -5))

net_list <- lapply(paste(path, net_names, sep = ""), read.table, header = TRUE, sep = "\t", check.names = FALSE)

sp_struct <- data.frame(ID = character(0), Species = character(0), MDim = numeric(0), Degree = numeric(0), GC_Aff = numeric(0), Fiedler = numeric(0)) # dataframe to be filled

# --------------- NODE AFFILIATION TO COMPONENTS AND SECTORS -------------- 

# position and ID of inverted networks
#invposit <- c(1,3,10,11,23,35,51,52,53,61,72,79,84,93,97,99,105,107,108,137,158,165,168,172,178,180,181,183,184,189,191,194,198,214,219,221,224,229,230,233,234,242,252,257,258,262,266,267,273,274,289,292,294,295,297,301,305,322,323,324,327,331,332,346,353,360,362,365,369,370); invID <- dataset$ID[invposit]

for (i in 1:length(net_names)) {

  ID <- dataset[i, "ID"]
  IntType <- dataset[i, "IntType"]
  M <- net_list[[i]]; M[M > 0] <- 1
  
  ID <- str_sub(net_names[i], end = -5)
  Species <- c(rownames(M), colnames(M))
  MDim <- c(rep(1, dim(M)[1]), rep(2, dim(M)[2]))
  Degree <- c(rowSums(M), colSums(M))
  fiedler <- component_sp(M) # returns the 
  #if(ID %in% invID) {fiedler$Fiedler <- -fiedler$Fiedler}
  
  # Creates matrices of networks' largest components for modularity analysis
  lcs <- unique(na.omit(fiedler$GC_aff)) # largest net components
  for (j in 1:length(lcs)) {
    lcsp <- fiedler[which(fiedler$GC_aff == j), "sp"]
    Mlc <- as.matrix(M[rownames(M) %in% lcsp, colnames(M) %in% lcsp])
    if (nrow(Mlc) != 1 & ncol(Mlc) != 1) {
      #write.table(Mlc, paste("./Modularity/Matrices/", "M_", dataset[i, "ID"], "_GC", j, ".txt", sep = ""), sep = "\t", row.names = FALSE, col.names = FALSE)
      SpList <- data.frame(Sp = c(rownames(Mlc), colnames(Mlc)), MDim = c(rep(1, nrow(Mlc)), rep(2, ncol(Mlc))))
      #write.table(SpList, paste("./Modularity/Sp_Names/", "Sp_", dataset[i, "ID"], "_GC", j, ".txt", sep = ""), sep = "\t", row.names = FALSE)
    }
  }
  
  df <- data.frame(ID = rep(ID, sum(dim(M)[1], dim(M)[2])), Species = Species, MDim = MDim, Degree = Degree, GC_Aff = fiedler$GC_aff, Fiedler = fiedler$Fiedler)
  sp_struct <- rbind(sp_struct, df)
}

head(sp_struct); dim(sp_struct)

# ---------------------- NODE AFFILIATION TO MODULES ----------------------

sp_struct$Mod_Aff <- NA

path_members <- "./Modularity/Matrices/SimulatedAnnealing/"
members <- grep("MEMBERS", dir(path_members), value = TRUE)

path_sp <- "./Modularity/Sp_Names/"
spnames <- dir(path_sp)

all(str_sub(members, start = 11) == str_sub(spnames, start = 4))

for (i in 1:length(spnames)) {
  
  ID <- str_sub(members[i], start = 11, end = -9)
  m <- read.table(paste(path_members, members[i], sep = ""), sep = "\t", header = TRUE)
  s <- read.table(paste(path_sp, spnames[i], sep = ""), sep = "\t", header = TRUE)
  mdim <- table(s$MDim)
  order <- c(paste("R", 1:mdim[1], sep = ""), paste("C", 1:mdim[2], sep = ""))
  m <- m[match(order, m$Node),]
  s$Mod <- m$Module
  sp_struct[intersect(which(sp_struct$ID == ID), which(sp_struct$Species %in% s$Sp)), "Mod_Aff"] <- s$Mod # CHECKED
  
}

head(sp_struct); dim(sp_struct)
#write.table(sp_struct, "./Outputs/01_Node-Level_Struct.txt", sep = "\t", row.names = FALSE)