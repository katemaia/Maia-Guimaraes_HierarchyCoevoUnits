# ------------------ Net & Node-level data: add Subgraph ------------------

# Author: Kate P Maia
# Checked: 12/2023

# Adding subgraph to network-level data: Computes the number of motif-6 within network's largest components, sectors and modules and tests motif prevalence with null model (Network-level motif analysis). Adds results to network-level data.
# Adding subgraph to node-level data: Counts the number of motifs each species takes part in. Uses the bmotif package.

# --------------------- Loading library, code and data --------------------

library(stringr)
library(bmotif)

dataset <- read.table("./Data/Project_Dataset.txt", header = T, sep = "\t", stringsAsFactors = F)

sp_struct <- read.table("./Outputs/01_Node-Level_Struct.txt", header = T, stringsAsFactors = F)
net_struct <- read.table("./Outputs/02_Net-Level_Struct.txt", header = T, stringsAsFactors = F)

dataset <- dataset[match(net_struct$ID, dataset$ID),]
all(dataset$ID == net_struct$ID)
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

# ----------------- Adding subgraph to network-level data -----------------

# ----- subgraph in module and sector -----

net_m6 <- data.frame(ID = str_sub(net_names, end = -5), n_m6 = NA, n_m6_mod = NA, mean_m6_mod = NA, sd_m6_mod = NA, p_m6_mod = NA, n_m6_sect = NA, mean_m6_sect = NA, sd_m6_sect = NA, p_m6_sect = NA)

nnull <- 100
for (i in 1:length(net_names)) {
  
  print(i); M <- net_list[[i]]; ID <- str_sub(net_names[i], end = -5)
  
  sp_M <- sp_struct[sp_struct$ID == ID, c("ID", "Species", "MDim", "Fiedler", "GC_Aff", "Mod_Aff")]
  sp_M <- sp_M[!is.na(sp_M$GC_Aff) & sp_M$GC_Aff == 1,] # sp in largest component
  sp_M$Mod_Aff <- sp_M$Mod_Aff + 1
  
  M <- as.matrix(M[rownames(M) %in% sp_M$Species, colnames(M) %in% sp_M$Species])
  if (any(c(rownames(M), colnames(M)) != sp_M$Species)) {print("ERROR")}
  
  # Motif count in the LARGEST COMPONENT
  nm6_M <- mcount(M, FALSE, FALSE, FALSE, FALSE)[6, "frequency"]
  
  # Motif count inside SECTORS
  sp_sect <- tapply(sp_M$Species, INDEX = sp_M$Fiedler, list) # sp in each sector
  M_sects <- lapply(sp_sect, function(x, M) M <- M[rownames(M) %in% x, colnames(M) %in% x], M) # M per sector
  M_sects <- lapply(M_sects, as.matrix)
  nm6_sects <- lapply(M_sects, mcount, FALSE, FALSE, FALSE, FALSE)
  nm6_sects <- lapply(nm6_sects, function(x) x[6,"frequency"]) # n m6 per sector
  
  # Motif count inside MODULES
  sp_mods <- tapply(sp_M$Species, INDEX = sp_M$Mod_Aff, list) # sp in each module
  M_mods <- lapply(sp_mods, function(x, M) M <- M[rownames(M) %in% x, colnames(M) %in% x], M) # M per module
  M_mods <- lapply(M_mods, as.matrix)
  nm6_mods <- lapply(M_mods, mcount, FALSE, FALSE, FALSE, FALSE)
  nm6_mods <- lapply(nm6_mods, function(x) x[6,"frequency"]) # n m6 per module
  
  # temporary result dataframes
  df_sects <- t(as.data.frame(unlist(nm6_sects)))
  df_mods <- t(as.data.frame(unlist(nm6_mods)))
  rownames(df_sects) <- rownames(df_mods) <- NULL
  df_sects <- cbind(df_sects, betw = (nm6_M - sum(df_sects)))
  df_mods <- cbind(df_mods, betw = (nm6_M - sum(df_mods)))
  
  for (j in 2:(nnull + 1)) {
    
    sp_M$Mod_Aff <- unlist(tapply(sp_M$Mod_Aff, sp_M$MDim, sample))
    sp_M$Fiedler <- unlist(tapply(sp_M$Fiedler, sp_M$MDim, sample))
    
    # Motif count inside SECTORS
    sp_sect <- tapply(sp_M$Species, INDEX = sp_M$Fiedler, list) # sp in each sector
    M_sects <- lapply(sp_sect, function(x, M) M <- M[rownames(M) %in% x, colnames(M) %in% x], M)
    M_sects <- lapply(M_sects, as.matrix)
    nm6_sects <- lapply(M_sects, mcount, FALSE, FALSE, FALSE, FALSE)
    nm6_sects <- lapply(nm6_sects, function(x) x[6,"frequency"]) # n m6 per sector
    
    # Motif count inside MODULES
    sp_mods <- tapply(sp_M$Species, INDEX = sp_M$Mod_Aff, list) # sp in each module
    M_mods <- lapply(sp_mods, function(x, M) M <- M[rownames(M) %in% x, colnames(M) %in% x], M)
    M_mods <- lapply(M_mods, as.matrix)
    nm6_mods <- lapply(M_mods, mcount, FALSE, FALSE, FALSE, FALSE)
    nm6_mods <- lapply(nm6_mods, function(x) x[6,"frequency"]) # n m6 per module
    
    # add resul to temp result data.frames
    df_s <- t(as.data.frame(unlist(nm6_sects)))
    df_m <- t(as.data.frame(unlist(nm6_mods)))
    rownames(df_s) <- rownames(df_m) <- NULL
    df_sects <- rbind(df_sects, cbind(df_s, betw = (nm6_M - sum(df_s))))
    df_mods <- rbind(df_mods, cbind(df_m, betw = (nm6_M - sum(df_m))))
  }
  
  result <- data.frame(m = rowSums(df_mods[,-ncol(df_mods)]), s = rowSums(df_sects[,-ncol(df_sects)]))
  
  net_m6[net_m6$ID == ID,"n_m6"] <- nm6_M
  net_m6[net_m6$ID == ID, c("n_m6_mod", "n_m6_sect")] <- result[1,] # empirical
  net_m6[net_m6$ID == ID, c("mean_m6_mod", "mean_m6_sect")] <- round(apply(result[-1,], 2, mean), 3)
  net_m6[net_m6$ID == ID, c("sd_m6_mod", "sd_m6_sect")] <- round(apply(result[-1,], 2, sd), 3)
  net_m6[net_m6$ID == ID, c("p_m6_mod", "p_m6_sect")] <- apply(result, 2, function(x) {sum(x[-1] >= x[1])/nnull})
}

net_struct$ID[!net_struct$ID %in% net_m6$ID] # missing: multLC and starLC
head(net_struct); dim(net_struct)

# adds net_m6 new info into net_struct
net_struct <- cbind(net_struct[,-18], net_m6[match(net_struct$ID, net_m6$ID), -1])

# ----- subgraph in module-sector (double hierarchy) -----

path <- "./Outputs/m6-Profiles/"; prf_names <- dir(path)
path <- "./Data/Matrices/"; net_names <- dir(path)
net_names <- net_names[net_names %in% gsub("prf_", "", prf_names)]
all(net_names == gsub("prf_", "", prf_names))

net_list <- lapply(paste("./Data/Matrices/", net_names, sep = ""), read.table, header = T, sep = "\t", check.names = F)
prf_list <- lapply(paste("./Outputs/m6-Profiles/", prf_names, sep = ""), read.table, header = T, sep = "\t")

net_m6 <- data.frame(ID = str_sub(net_names, end = -5), n_m6_ms = NA, mean_m6_ms = NA, sd_m6_ms = NA, p_m6_ms = NA)

nnull <- 100
for (i in 1:length(net_names)) {
  
  print(i)
  ID <- str_sub(net_names[i], end = -5); M <- net_list[[i]]
  
  sp_M <- sp_struct[sp_struct$ID == ID, c("ID", "Species", "MDim", "Fiedler", "GC_Aff", "Mod_Aff")]
  sp_M <- sp_M[!is.na(sp_M$GC_Aff) & sp_M$GC_Aff == 1,]
  sp_M$Mod_Aff <- sp_M$Mod_Aff + 1
  
  if (any(rownames(M) %in% colnames(M)) | any(colnames(M) %in% rownames(M))) {print("ERROR")}
  
  M <- M[rownames(M) %in% sp_M$Species, colnames(M) %in% sp_M$Species]
  if (any(c(rownames(M), colnames(M)) != sp_M$Species)) {print("ERROR")}
  
  profile <- prf_list[[i]]
  
  if (nrow(profile) == 0) {net_m6[i, 2:5] <- c(0, 0, 0, 1)} else {
    
    prof_mod <- prof_sect <- profile
    
    # network-level motif analysis
    prof_mod$R1 <- sp_M[match(profile$R1, sp_M$Species), "Mod_Aff"]
    prof_mod$R2 <- sp_M[match(profile$R2, sp_M$Species), "Mod_Aff"]
    prof_mod$C1 <- sp_M[match(profile$C1, sp_M$Species), "Mod_Aff"]
    prof_mod$C2 <- sp_M[match(profile$C2, sp_M$Species), "Mod_Aff"]
    
    prof_sect$R1 <- sp_M[match(profile$R1, sp_M$Species), "Fiedler"]
    prof_sect$R2 <- sp_M[match(profile$R2, sp_M$Species), "Fiedler"]
    prof_sect$C1 <- sp_M[match(profile$C1, sp_M$Species), "Fiedler"]
    prof_sect$C2 <- sp_M[match(profile$C2, sp_M$Species), "Fiedler"]
    
    # null model: shuffles module/sector affiliation
    result <- rep(NA, (nnull + 1))
    motinmod <- which(apply(apply(prof_mod, 1, diff) == 0, 2, all)) # all mot sp in one mod
    motinsect <- which(apply(apply(prof_sect, 1, diff) == 0, 2, all)) # all mot sp in one sect
    result[1] <- length(intersect(motinmod, motinsect))
    
    sp_M$mstag <- paste(sp_M$Mod_Aff, sp_M$Fiedler, sep = "_")
    for (j in 2:(nnull + 1)) {
      
      sp_M$mstag <- unlist(tapply(sp_M$mstag, sp_M$MDim, sample))
      sp_M$Mod_Aff <- as.numeric(unlist(lapply(str_split(sp_M$mstag, "_"), function(x) x[1])))
      sp_M$Fiedler <- as.numeric(unlist(lapply(str_split(sp_M$mstag, "_"), function(x) x[2])))
      
      prof_mod <- prof_sect <- profile
      prof_mod$R1 <- sp_M[match(profile$R1, sp_M$Species), "Mod_Aff"]
      prof_mod$R2 <- sp_M[match(profile$R2, sp_M$Species), "Mod_Aff"]
      prof_mod$C1 <- sp_M[match(profile$C1, sp_M$Species), "Mod_Aff"]
      prof_mod$C2 <- sp_M[match(profile$C2, sp_M$Species), "Mod_Aff"]
      
      prof_sect$R1 <- sp_M[match(profile$R1, sp_M$Species), "Fiedler"]
      prof_sect$R2 <- sp_M[match(profile$R2, sp_M$Species), "Fiedler"]
      prof_sect$C1 <- sp_M[match(profile$C1, sp_M$Species), "Fiedler"]
      prof_sect$C2 <- sp_M[match(profile$C2, sp_M$Species), "Fiedler"]
      
      motinmod <- which(apply(apply(prof_mod, 1, diff) == 0, 2, all)) # prof_mod[motinmod,]
      motinsect <- which(apply(apply(prof_sect, 1, diff) == 0, 2, all)) # prof_sect[motinsect,]
      result[j] <- length(intersect(motinmod, motinsect))
    }
    
    net_m6[net_m6$ID == ID, 2:4] <- c(result[1], mean(result[-1]), sd(result[-1]))
    net_m6[net_m6$ID == ID, "p_m6_ms"] <- sum(result[-1] >= result[1])/nnull
  }
}

head(net_m6); dim(net_m6)

net_struct <- cbind(net_struct, net_m6[match(net_struct$ID, net_m6$ID), -1])
head(net_struct); dim(net_struct)

#write.table(net_struct, "./Outputs/04_Net-Level_Struct.txt", sep = "\t", row.names = FALSE)

# ------------------ Adding subgraph to node-level data ------------------

# Takes a long time as it uses a jack-knifing procedure.

sp_struct$m6 <- NA
for (i in 1:length(net_names)) {
  
  print(i); M <- net_list[[i]]; ID <- str_sub(net_names[i], end = -5)

  sp_M <- sp_struct[sp_struct$ID == ID, c("ID", "Species", "MDim", "Fiedler", "GC_Aff", "Mod_Aff")]
  sp_M <- sp_M[!is.na(sp_M$GC_Aff) & sp_M$GC_Aff == 1,]
  sp_M$Mod_Aff <- sp_M$Mod_Aff + 1
  
  M <- as.matrix(M[rownames(M) %in% sp_M$Species, colnames(M) %in% sp_M$Species])
  if (any(c(rownames(M), colnames(M)) != sp_M$Species)) {print("ERROR")}
  
  # motif count in complete network (largest component)
  nm6_M <- mcount(M, FALSE, FALSE, FALSE, FALSE)[6, "frequency"]
  
  sp_nm6 <- c()
  if (nrow(M) > 2) {
    for (j in 1:nrow(M)) {sp_nm6 <- c(sp_nm6, nm6_M - mcount(M[-j,], F, F, F, F)[6,"frequency"])}
  } else {sp_nm6 <- c(sp_nm6, rep(nm6_M, nrow(M)))} # if nrow <= 2, nrow sp will be in nm6_M
  if (ncol(M) > 2) {
    for (j in 1:ncol(M)) {sp_nm6 <- c(sp_nm6, nm6_M - mcount(M[,-j], F, F, F, F)[6,"frequency"])}
  } else {sp_nm6 <- c(sp_nm6, rep(nm6_M, ncol(M)))} # if ncol <= 2, ncol sp will be in nm6_M
  
  if (any(sp_struct[!is.na(sp_struct$GC_Aff) & sp_struct$GC_Aff == 1 & sp_struct$ID == ID, "Species"] != 
    c(rownames(M), colnames(M)))) {print("ERROR")}
  sp_struct[!is.na(sp_struct$GC_Aff) & sp_struct$GC_Aff == 1 & sp_struct$ID == ID, "m6"] <- sp_nm6
}

sp_struct <- sp_struct[, -8] # removes interaction type
#write.table(sp_struct, "./Outputs/04_Node-Level_Struct.txt", sep = "\t", row.names = FALSE)