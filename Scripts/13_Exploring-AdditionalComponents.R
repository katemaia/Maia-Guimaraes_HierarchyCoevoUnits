# --------------- Exploring networks' additional components ---------------

# Author: Kate P Maia
# Checked: 08/2024

# Explores networks' additional components, those which are not the largest. Code is divided in:
# Loading library, code and data
# Creating component level data on mult-comp networks: creates data on each components of the 143 networks with multiple components. Data includes: network identification (ID), component number (CNo), component size (CSize), whether that is network's largest component (isLC, 1 if it is and 0 if not), component's number species in each set and of interactions (NRow, NCol, NInt).
# Exploring and plotting component level data on mult-comp networks

# --------------------- Loading library, code and data --------------------

library(tidyverse)
library(igraph)
library(cowplot)
library(scales)
source("./Scripts/Functions/inc2adj.R")

net_struct <- read.table("./Outputs/05_Net-Level_Struct.txt", header = TRUE, sep = "\t")

multComp <- net_struct[net_struct$NComp != 1, 1]; length(multComp) # multi-component networks (N = 143)
net_list <- lapply(paste("./Data/Matrices/", multComp, ".txt", sep = ""), read.table, header = TRUE, sep = "\t", check.names = FALSE); length(net_list)

# ---------- Creating component level data on mult-comp networks ----------

resultnet <-  data.frame(ID = character(0), CSize = numeric(0), isLC = numeric(0), NRow = numeric(0), NCol = numeric(0), NInt = numeric(0))

for (i in 1:length(multComp)) { # for each of the mult-component networks
  
  M <- net_list[[i]]; M[M > 0] <- 1 # binary
  ID <- multComp[i]; print(ID)
  
  A <- inc2adj(M) # adjacency matrix
  g <- graph_from_adjacency_matrix(A, "undirected")
  clu <- igraph::components(g)
  if (any(clu$Species != c(rownames(M), colnames(M)))) {print("ERROR: species order")}
  
  df <- data.frame(NRow = numeric(0), NCol = numeric(0), NInt = numeric(0))
  for (j in 1:clu$no){ # structural description of each component
    memb <- clu$membership[clu$membership == j]
    MC <- as.matrix(M[rownames(M) %in% names(memb), colnames(M) %in% names(memb)])
    df <- rbind(df, data.frame(NRow = nrow(MC), NCol = ncol(MC), NInt = sum(MC)))
  }
  
  net <-  data.frame(ID, CNo = 1:clu$no, CSize = clu$csize, isLC = ifelse(clu$csize == max(clu$csize), 1, 0))
  net <-  cbind(net, df) # add info on component structure
  resultnet <- rbind(resultnet, net)
  
}

all(resultnet$CSize == resultnet$NRow + resultnet$NCol)

# ------ Exploring and plotting component data on mult-comp networks ------ 

net_struct$IntType <- factor(net_struct$IntType, levels = c("Host-Parasite", "Phage-Bacteria", "Plant-Herbivore", "Food Webs", "Pollination", "Seed Dispersal", "Anemone-Fish", "Plant-Ant"))
net_struct$Code <- factor(net_struct$IntType)
levels(net_struct$Code) <- c("Pa-Ho", "Ph-B", "Her-P", "Pr-Pr", "Pol-P", "SD-P", "F-Ane", "Ant-P")

pal <- hue_pal()(8) # palette

# Single vs Mult Component networks per interaction type
pa <- net_struct %>% select(ID, Code, NComp) %>% mutate(NComp = NComp == 1) %>%
  mutate(NComp = factor(NComp, levels = c("TRUE", "FALSE"))) %>% 
  group_by(Code, NComp) %>% count() %>% 
  ggplot(aes(x = NComp, y = n, fill = Code)) +
  geom_bar(position = "stack", stat = "identity") + 
  xlab("Number of components") + ylab("Number of networks") + 
  scale_fill_manual(values = pal) + scale_x_discrete(labels = c("One", "More")) + 
  theme_cowplot(); pa

# Checking if net_struct matches with resultnet
all(table(resultnet$ID) == net_struct[net_struct$ID %in% multComp, "NComp"]) # no of components
all(unique(resultnet[resultnet$isLC == 1, c("ID", "CSize")])[,2] == net_struct[net_struct$ID %in% multComp, "LargeCompSize"]) # size of largest component

# Understanding size of networks' additional components
resultnet <- as_tibble(resultnet); dim(resultnet) # 778 components in 143 networks

starLC <- c("IZ_PA_DIM10", "IZ_PA_PA10", "M_AF_002_07", "M_AF_002_16", "M_PL_061_33")
multLC <- c("Herb_CuevasReyes", "IZ_PA_DIM100", "IZ_PA_PAcamp", "M_AF_002_02", "M_AF_002_11")
all(c(starLC, multLC) %in% multComp) # including 10 networks excluded from analysis
resultnet %>% filter(ID %in% multLC) %>% print(n = 46) # 46 from multLC, 2-5 species
resultnet %>% filter(ID %in% starLC) %>% print(n = 21) # 21 from starLC, 2-6 species

resultnet %>% filter(isLC != 1) # 601/778 are additional (non-largest) components
resultnet %>% filter(isLC != 1) %>% filter(CSize == 2) # 394 (66%) have 2 species
resultnet %>% filter(isLC != 1) %>% filter(CSize == 3) # 99 (16%) have 3 species
resultnet %>% filter(isLC != 1) %>% filter(CSize >= 4) %>%
  filter(NRow == 1 | NCol == 1) # 54 (9%) are star components with 4+ species

# Plotting additional component results (601 components)
addcomp <- resultnet %>% filter(isLC != 1) %>% mutate(CSizeFct = NA)
addcomp[addcomp$CSize == 2, "CSizeFct"] <- 2
addcomp[addcomp$CSize == 3, "CSizeFct"] <- 3
addcomp[addcomp$CSize > 3 & addcomp$NRow == 1, "CSizeFct"] <- 5 # 5 is cue for star 
addcomp[addcomp$CSize > 3 & addcomp$NCol == 1, "CSizeFct"] <- 5 # 5 is cue for star
addcomp[is.na(addcomp$CSizeFct), "CSizeFct"] <- 7 # 7 is cue for other
addcomp$CSizeFct <- factor(addcomp$CSizeFct, levels = c(2, 3, 5, 7))
levels(addcomp$CSizeFct) <- c("2 Sp", "3 Sp", "Star", "Other")
addcomp <- left_join(addcomp, net_struct[,c(1, 3, 4, 35)], by = "ID")

pb <- addcomp %>% group_by(Code, CSizeFct) %>% count() %>% 
  ggplot(aes(x = CSizeFct, y = n, fill = Code)) +
  geom_bar(position = "stack", stat = "identity") + 
  xlab("") + ylab("Number of additional \n components (not largest)") + 
  scale_fill_manual(values = pal[-4]) + theme_cowplot() +
  theme(legend.position = "none", axis.title.y = element_text(size = 12)); pb
pb <- pb + annotate("text", x = 1:4, y = c(394, 99, 54, 54) + 12, label = c("66%", "16%", "9%", "9%"))

resultnet %>% filter(isLC != 1) %>% filter(CSize >= 4) %>%
  filter(NRow != 1 & NCol != 1) # 54 are Other (4+ species and not stars)
resultnet %>% filter(isLC != 1) %>% filter(CSize >= 4) %>% 
  filter(NRow != 1 & NCol != 1) %>% filter(CSize == 4) # 15 2x2, 4 motifs

# Plotting "Other" component results (54 components)
othercomp <- addcomp %>% filter(CSizeFct == "Other") 
othercomp <- othercomp %>% mutate(CSizeP = CSize/(NRow.y + NCol.y)) # size as proportion of network species

pc <- othercomp %>% ggplot(aes(x = CSize, y = CSizeP, color = Code)) + 
  geom_point(size = 3, alpha = 0.7) + 
  scale_color_manual(values = pal[c(1, 3, 5, 8)]) + 
  xlab("Component size (species number)") + 
  ylab("Component size\n (network species proportion)") + 
  scale_x_continuous(breaks = c(4, 8, 12, 16), limits = c(3,16)) + 
  theme_cowplot() + theme(legend.position = "none", axis.title.x = element_text(size = 11.5), axis.title.y = element_text(size = 12)); pc

plot_grid(pa, pb, pc, labels = c("a)", "b)", "c)"), nrow = 1, label_size = 12, label_fontfamily = "", label_fontface = "plain")
