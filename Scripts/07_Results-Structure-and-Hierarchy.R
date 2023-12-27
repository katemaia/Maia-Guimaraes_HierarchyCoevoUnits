# --------------- Manuscript figures: groups and hierarchy ----------------

# Author: Kate P Maia
# Checked: 12/2023

# Reads data generated in previous scripts to generate manuscript results on group structure of ecological networks and hierarchical structure. 

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

net_struct %>% filter(!is.na(n_m6)) %>% # Number of subgraphs
  group_by(IntType) %>% summarise(mean = mean(n_m6), sd = sd(n_m6), min = min(n_m6), max = max(n_m6))

net_struct %>% filter(!is.na(n_m6)) %>% 
  group_by(IntSign) %>% summarise(mean = mean(n_m6), sd = sd(n_m6))
anova(lm(n_m6 ~ IntSign, data = net_struct))

# -------------------------------------------------------------------------
# ------------------ Results on module-sector congruence ------------------

net_struct %>% group_by(IntType) %>% # Module-sector congruence
  summarise(mean = mean(MS_congr), sd = sd(MS_congr), min = min(MS_congr), max = max(MS_congr))

# -------------------------------------------------------------------------
# ------------ Results on subgraph and module/sector congruence -----------

# average number of subgraphs in sectors across networks
net_struct %>% select(Code, ID, n_m6, n_m6_sect) %>% # all network types (N = 366)
  mutate(In = n_m6_sect/n_m6) %>% filter(!is.nan(In)) %>% # 345 remaining as NaNs are 0/0 (21 had none)
  summarise(M = mean(In)*100, SD = sd(In)*100)

# proportion of networks with a significant number of subgraphs in sectors (null model)
net_struct %>% select(Code, p_m6_sect) %>%
  mutate(p_m6_sect = ifelse(p_m6_sect <= 0.05, TRUE, FALSE)) %>% 
  summarise(M = mean(p_m6_sect)*100)

# average number of subgraphs in modules across networks
net_struct %>% select(Code, ID, n_m6, n_m6_mod) %>% # all network types (N = 366)
  mutate(In = n_m6_mod/n_m6) %>% filter(!is.nan(In)) %>% # 345 remaining as NaNs are 0/0 (21 had none)
  summarise(M = mean(In)*100, SD = sd(In)*100)

# proportion of networks with a significant number of subgraphs in modules (null model)
net_struct %>% select(Code, p_m6_mod) %>%
  mutate(p_m6_mod = ifelse(p_m6_mod <= 0.05, TRUE, FALSE)) %>% 
  summarise(M = mean(p_m6_mod)*100)

#*also excluded M_PL_62 as the null model uses profile instead of bmotif (need motif identity)
net_struct %>% select(Code, ID, n_m6, n_m6_ms) %>% # all network types (N = 366)
  mutate(In = n_m6_ms/n_m6) %>% filter(!is.na(In)) %>% # 344 remaining as NaNs are 0/0 (21 had none, + M_PL_62)
  summarise(M = mean(In)*100, SD = sd(In)*100)

# *also excluded M_PL_62 as the null model uses profile instead of bmotif (need motif identity)
net_struct %>% select(Code, p_m6_ms) %>%
  mutate(p_m6_ms = ifelse(p_m6_ms <= 0.05, TRUE, FALSE)) %>% filter(!is.na(p_m6_ms)) %>% 
  summarise(M = mean(p_m6_ms)*100)