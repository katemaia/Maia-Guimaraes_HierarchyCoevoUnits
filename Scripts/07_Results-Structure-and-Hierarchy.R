# --------------- Manuscript results: groups and hierarchy ----------------

# Author: Kate P Maia
# Checked: 12/2023

# PART1: Reads data generated in previous scripts to generate manuscript results on group structure of ecological networks and hierarchical structure. 
# PART2: Complements Table 1 (network structure - groups across scales) with z-scores to control for connectance and species richness. 

# --------------------- Loading library, code and data --------------------

library(tidyverse)

dataset <- read.table("./Data/Project_Dataset.txt", header = T, sep = "\t", stringsAsFactors = F)

net_struct <- read.table("./Outputs/05_Net-Level_Struct.txt", header = T, stringsAsFactors = F)
sp_struct <- read.table("./Outputs/04_Node-Level_Struct.txt", header = T, stringsAsFactors = F)

dataset <- dataset[match(net_struct$ID, dataset$ID),]; all(dataset$ID == net_struct$ID)
sp_struct <- sp_struct[order(sp_struct$ID),]; all(dataset$ID == levels(sp_struct$ID))

net_struct$IntType <- factor(net_struct$IntType, levels = c("Host-Parasite", "Phage-Bacteria", "Plant-Herbivore", "Food Webs", "Pollination", "Seed Dispersal", "Anemone-Fish", "Plant-Ant"))
net_struct$IntSign <- factor(dataset$IntSign)
net_struct$Code <- factor(net_struct$IntType)
levels(net_struct$IntType) <- c("Parasite-Host", "Phage-Bacteria", "Herbivore-Plant", "Predator-Prey", "Pollinator-Plant", "Seed Disperser-Plant", "Fish-Anemone", "Ant-Plant")
levels(net_struct$Code) <- c("Pa-Ho", "Ph-B", "Her-P", "Pr-Pr", "Pol-P", "SD-P", "F-Ane", "Ant-P")

sp_struct$IntType <- net_struct[match(sp_struct$ID, net_struct$ID), "IntType"]
sp_struct$IntSign <- net_struct[match(sp_struct$ID, net_struct$ID), "IntSign"]

dim(net_struct); length(unique(sp_struct$ID))

# -------------------------------------------------------------------------
# ---------------- Results on basic network-level structure ---------------

table(net_struct$IntSign); table(net_struct$IntType) # Table 1 (Sign, Type, N)

net_struct$S <- net_struct$NRow + net_struct$NCol # Network species richness
anova(lm(S ~ IntSign, data = net_struct))
aggregate(S ~ IntSign, data = net_struct, mean)
aggregate(S ~ IntSign, data = net_struct, sd)

anova(lm(Connectance ~ IntSign, data = net_struct)) # Network connectance
aggregate(Connectance ~ IntSign, data = net_struct, mean)
aggregate(Connectance ~ IntSign, data = net_struct, sd)

m_S <- aggregate(S ~ IntType, data = net_struct, mean)
sd_S <- aggregate(S ~ IntType, data = net_struct, sd)
m_C <- aggregate(Connectance ~ IntType, data = net_struct, mean)
sd_C <- aggregate(Connectance ~ IntType, data = net_struct, sd)
all(m_S$IntType == sd_S$IntType); all(m_S$IntType == m_C$IntType)
summary_df <- data.frame(Type = m_S[,1], m_S = m_S[,2], sdS = sd_S[,2], C = m_C[,2], sdC = sd_C[,2])
summary_df[,-1] <- round(summary_df[,-1], digits = 2); summary_df # Table 1 (S, Connectance)

# -------------------------------------------------------------------------
# --------------------- Results on component structure --------------------

sum(net_struct$NComp == 1); sum(net_struct$NComp == 1)/nrow(net_struct) # Single component networks

net_struct %>% mutate(NComp = NComp == 1) %>% 
  group_by(IntType) %>% count(NComp) %>% 
  pivot_wider(names_from = NComp, values_from = n) %>% 
  transmute(SingComp = `TRUE`, All = `FALSE` + `TRUE`)

anova(lm(NComp ~ IntSign, data = net_struct)) # Network fragmentation (number of components)
aggregate(NComp ~ IntSign, data = net_struct, mean)
aggregate(NComp ~ IntSign, data = net_struct, sd)

m_NComp <- aggregate(NComp ~ IntType, data = net_struct, mean)
sd_NComp <- aggregate(NComp ~ IntType, data = net_struct, sd)
all(m_NComp$IntType == sd_NComp$IntType)
m_NComp$NComp <- round(m_NComp$NComp, digits = 2)
cbind(m_NComp, sd_NComp = round(sd_NComp$NComp, digits = 2)) # Table 1 (Number of components)

net_struct$SizeLG <- round(net_struct$LargeCompSize/net_struct$S, digits = 2) # Largest component species richness
anova(lm(SizeLG ~ IntSign, data = net_struct))
aggregate(SizeLG ~ IntSign, data = net_struct, mean)
aggregate(SizeLG ~ IntSign, data = net_struct, sd)

m_SizeLG <- aggregate(SizeLG ~ IntType, data = net_struct, mean)
sd_SizeLG <- aggregate(SizeLG ~ IntType, data = net_struct, sd)
all(m_SizeLG$IntType == sd_SizeLG$IntType)
m_SizeLG$SizeLG <- round(m_SizeLG$SizeLG, digits = 2)
cbind(m_SizeLG, sd_SizeLG = round(sd_SizeLG$SizeLG, digits = 2)) # Table 1 (Largest component)

# -------------------------------------------------------------------------
# ----------------- Results on modules (largest component) ----------------

# REMOVING multLC & starLC for further analysis (within larget components)
multLC <- c("Herb_CuevasReyes", "IZ_PA_DIM100", "IZ_PA_PAcamp", "M_AF_002_02", "M_AF_002_11")
starLC <- c("IZ_PA_DIM10", "IZ_PA_PA10", "M_AF_002_07", "M_AF_002_16", "M_PL_061_33")

net_struct <- net_struct[!net_struct$ID %in% c(multLC, starLC),]; dim(net_struct) # N = 366

anova(lm(SA_Modularity ~ IntSign, data = net_struct)) # Modularity
net_struct %>% select(IntSign, SA_Modularity) %>% group_by(IntSign) %>% 
  summarise(mean = mean(SA_Modularity), sd = sd(SA_Modularity))

anova(lm(SA_NMod ~ IntSign, data = net_struct)) # Number of modules
net_struct %>% select(IntSign, SA_NMod) %>% group_by(IntSign) %>% 
  summarise(mean = mean(SA_NMod), sd = sd(SA_NMod))

net_struct %>% group_by(IntType) %>% 
  summarise(mean = mean(SA_NMod), sd = sd(SA_NMod)) # Table 1 (Number of modules)

# -------------------------------------------------------------------------
# ---------------- Results on subgraphs (largest component) ---------------

net_struct %>% filter(is.na(n_m6))
net_struct %>% filter(n_m6 == 0) %>% group_by(IntType) %>% count() # Networks without subgraphs

net_struct %>% group_by(IntType) %>% # Number of subgraphs
  summarise(mean = mean(n_m6), sd = sd(n_m6), min = min(n_m6), max = max(n_m6))

net_struct %>% group_by(IntSign) %>% 
  summarise(mean = mean(n_m6), sd = sd(n_m6))
anova(lm(n_m6 ~ IntSign, data = net_struct))

# -------------------------------------------------------------------------
# ------------------ Results on module-sector congruence ------------------

net_struct %>% group_by(IntType) %>% # Module-sector congruence
  summarise(mean = mean(MS_congr), sd = sd(MS_congr), min = min(MS_congr), max = max(MS_congr))

net_struct %>% mutate(MS_pval = MS_pval <= 0.05) %>% 
  group_by(IntType, MS_pval) %>% count() %>% arrange(MS_pval)

# -------------------------------------------------------------------------
# ------------ Results on subgraph and module/sector congruence -----------

# average number of subgraphs in sectors across networks
net_struct %>% select(Code, ID, n_m6, n_m6_sect) %>% # all network types (N = 366)
  mutate(In = n_m6_sect/n_m6) %>% filter(n_m6 != 0) %>% # 345 remaining as NaNs are 0/0 (21 had none)
  summarise(M = mean(In)*100, SD = sd(In)*100)

# proportion of networks with a significant number of subgraphs in sectors (null model)
net_struct %>% filter(n_m6 != 0) %>% # 345 remaining as NaNs are 0/0 (21 had none)
  select(Code, p_m6_sect) %>% mutate(p_m6_sect = ifelse(p_m6_sect <= 0.05, TRUE, FALSE)) %>% 
  summarise(M = mean(p_m6_sect)*100)

# average number of subgraphs in modules across networks
net_struct %>% select(Code, ID, n_m6, n_m6_mod) %>% # all network types (N = 366)
  mutate(In = n_m6_mod/n_m6) %>% filter(n_m6 != 0) %>% # 345 remaining as NaNs are 0/0 (21 had none)
  summarise(M = mean(In)*100, SD = sd(In)*100)

# proportion of networks with a significant number of subgraphs in modules (null model)
net_struct %>% filter(n_m6 != 0) %>% # 345 remaining as NaNs are 0/0 (21 had none)
  select(Code, p_m6_mod) %>% mutate(p_m6_mod = ifelse(p_m6_mod <= 0.05, TRUE, FALSE)) %>% 
  summarise(M = mean(p_m6_mod)*100)

#*also excluded M_PL_62 as the null model uses profile instead of bmotif (need motif identity)
net_struct %>% select(Code, ID, n_m6, n_m6_ms) %>% # all network types (N = 366)
  mutate(In = n_m6_ms/n_m6) %>% filter(n_m6 != 0) %>% filter(ID != "M_PL_062") %>% # 344 remaining as NaNs are 0/0 (21 had none, + M_PL_62)
  summarise(M = mean(In)*100, SD = sd(In)*100)

# *also excluded M_PL_62 as the null model uses profile instead of bmotif (need motif identity)
net_struct %>% filter(n_m6 != 0) %>% # 345 remaining as NaNs are 0/0 (21 had none) 
  select(Code, p_m6_ms) %>%
  mutate(p_m6_ms = ifelse(p_m6_ms <= 0.05, TRUE, FALSE)) %>% filter(!is.na(p_m6_ms)) %>% 
  summarise(M = mean(p_m6_ms)*100)

# -------------- PART 2: Complementing Table 1 with z-scores --------------

# Author: Kate P Maia
# Checked: 04/2024

# --------------------- Loading library, code and data --------------------

library(tidyverse)
library(bipartite)
library(igraph)
library(bmotif)
#source(./functions/inc2adj.R)

net_struct <- read.table("./Outputs/05_Net-Level_Struct.txt", header = T, sep = "\t")

path <- "./Data/Matrices/"; net_names <- dir(path)
net_list <- lapply(paste(path, net_names, sep = ""), read.table, header = TRUE, sep = "\t", check.names = FALSE)

path <- "./Modularity/Matrices/"; gc_names <- grep("M_", dir(path), value = TRUE)
# REMOVING multLC & starLC (multiple largest components, star-shaped largest components)
multLC <- c("Herb_CuevasReyes", "IZ_PA_DIM100", "IZ_PA_PAcamp", "M_AF_002_02", "M_AF_002_11")
starLC <- c("IZ_PA_DIM10", "IZ_PA_PA10", "M_AF_002_07", "M_AF_002_16", "M_PL_061_33")
gc_names <- gc_names[-grep(paste(c(starLC, multLC), collapse = "|"), gc_names)]
gc_list <- lapply(paste(path, gc_names, sep = ""), read.table, header = FALSE, sep = "\t", check.names = FALSE)

nrow(net_struct); length(net_names); length(net_list) # 376 networks
length(gc_names); length(gc_list) # 366 networks largest components

# ------------------- computing components and subgraphs ------------------

table1 <- data.frame(ID = character(0), NC = numeric(0), MNC = numeric(0), SDNC = numeric(0), Nm6 = numeric(0), MNm6 = numeric(0), SDNm6 = numeric(0))

for (i in 1:length(net_names)) { # for each network 
  
  print(i)
  ID <- gsub(".txt", "", net_names[i])
  M <- as.matrix(net_list[[i]]); M[M != 0] <- 1 # binary incidence
  flag <- ID %in% c(starLC, multLC, "Herb_Pearse&Alterm", "M_PL_062") # eliminated from subgraph analysis
  
  # ----- Network components ----- #
  
  Mnull <- nullmodel(M, N = 100, method = "shuffle.web") # null model
  conn <- any(unlist(lapply(Mnull, sum)) != sum(M))
  nrow <- any(unlist(lapply(Mnull, nrow)) != nrow(M))
  ncol <- any(unlist(lapply(Mnull, ncol)) != ncol(M))
  if (any(conn, nrow, ncol)) {print("ERROR IN NULL: M")} # confirms null model
  
  A <- inc2adj(M); Anull <- lapply(Mnull, inc2adj) # creates adjacency
  Anullg <- lapply(Anull, graph_from_adjacency_matrix, "undirected") # creates graphs
  
  NC <- count_components(graph_from_adjacency_matrix(A, "undirected")) # counting components
  NCn <- unlist(lapply(Anullg, count_components))
  compdf <- data.frame(NC, MNC = mean(NCn), SDNC = sd(NCn))
  
  # ----- Subgraphs in largest component ----- #
  
  if (flag) {m6df <- data.frame(Nm6 = NA, MNm6 = NA, SDNm6 = NA)
  } else {
    
    GCM <- as.matrix(gc_list[[grep(ID, gc_names)]]) # matrix of largest component
    Nm6 <- mcount(GCM, F, F, F, F)[6,] # counting subgraphs
    
    GCnull <- nullmodel(GCM, N = 100, method = "shuffle.web")
    conn <- any(unlist(lapply(GCnull, sum)) != sum(GCM))
    nrow <- any(unlist(lapply(GCnull, nrow)) != nrow(GCM))
    ncol <- any(unlist(lapply(GCnull, ncol)) != ncol(GCM))
    if (any(conn, nrow, ncol)) {print("ERROR IN NULL: GC")}
    
    Nm6n <- lapply(GCnull, mcount, F, F, F, F) # counting subgraphs in null model
    Nm6n <- do.call(rbind, lapply(Nm6n, function(x) x[6,]))
    m6df <- data.frame(Nm6 = Nm6$frequency, MNm6 = mean(Nm6n$frequency), SDNm6 = sd(Nm6n$frequency))
    
    # Saves null matrices to compute number of modules: Newman metric, Fast Grid algorithm, MODULAR PROGRAMME 
    nullnames <- paste0("./Outputs/null_GC_M/M_GC_", ID, "_N", 1:100, ".txt")
    #mapply(write.table, GCnull, nullnames, sep = "\t", row.names = FALSE, col.names = FALSE)
  }
  
  idf <- cbind(ID = ID, compdf, m6df); table1 <- rbind(table1, idf)
  
}

#write.table(table1, "./Outputs/07_No-of-Groups_NullModel.txt", sep = "\t")

# ------------------- computing zscore values for Table1 ------------------

table1 <- read.table("./Outputs/07_No-of-Groups_NullModel.txt", header = TRUE, sep = "\t")
modresults <- read.table("./Data/Null_GC_M/ResultsFG/OUT_MOD.txt", header = TRUE, sep = "\t") # number of modules: Newman metric, Fast Grid algorithm, MODULAR PROGRAMME 

# fixing network ID
modresults <- modresults %>% mutate(File = str_sub(File, start = 6, end = -5)) %>% 
  mutate(File = gsub(paste(paste0("_N", 1:100), collapse = "|"), "", File)) %>% 
  select(File, N.modules) %>% rename(ID = File, NMod = N.modules) %>% 
  group_by(ID) %>% summarise(MNM = mean(NMod), SDNM = sd(NMod))

# adding modularity results to table1
table1 <- left_join(table1, modresults, by = "ID") %>% select(1:4, MNM, SDNM, 5:7)

# adding IntType and number of modules to table1
table1 <- left_join(table1, net_struct[,c("ID", "IntType", "FG_NMod")], by = "ID") %>% select(1, 10, 2:4, 11, 5:9)

# calculating zscores for Table 1 of the manuscript
table1 <- table1 %>% mutate(ZComp = (NC - MNC)/SDNC, ZMod = (FG_NMod - MNM)/SDNM, ZSubg = (Nm6 - MNm6)/SDNm6)

# fixing interaction types
table1 <- table1 %>% mutate(IntType = factor(IntType, levels = c("Host-Parasite", "Phage-Bacteria", "Plant-Herbivore", "Food Webs", "Pollination", "Seed Dispersal", "Anemone-Fish", "Plant-Ant"))) %>%
  mutate(IntType = recode_factor(IntType, "Host-Parasite" = "Parasite-Host", "Phage-Bacteria" = "Phage-Bacteria", "Plant-Herbivore" = "Herbivore-Plant", "Food Webs" = "Predator-Prey", "Pollination" = "Pollinator-Plant", "Seed Dispersal" = "Seed Disperser-Plant", "Anemone-Fish" = "Fish-Anemone", "Plant-Ant" = "Ant-Plant")) %>% 
  mutate(SDNC = SDNC != 0, SDNM = SDNM != 0, SDNm6 = SDNm6 != 0)

# Number of Components  
table1 %>% filter(!SDNC) # 131 with sd = 0
table1 %>% filter(SDNC) %>% # 245 = 376 - 131
  group_by(IntType) %>% summarise(N = sum(SDNC), Mean = mean(ZComp), SD = sd(ZComp))

# Number of Modules  
table1 %>% filter(!SDNM) # 31 with sd = 0
table1 %>% filter(SDNM) %>% # 335 = 376 - 10 - 31 (10 starLC and multLC)
  group_by(IntType) %>% summarise(N = sum(SDNM), Mean = mean(ZMod), SD = sd(ZMod))

# Number of Subgraphs  
table1 %>% filter(!SDNm6) # 30 with sd = 0
table1 %>% filter(SDNm6) %>% # 336 = 376 - 10 - 30 (10 starLC and multLC)
  group_by(IntType) %>% summarise(N = sum(SDNm6), Mean = mean(ZSubg), SD = sd(ZSubg))
