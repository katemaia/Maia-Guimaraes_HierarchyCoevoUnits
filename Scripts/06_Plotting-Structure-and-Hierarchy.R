# --------------- Manuscript figures: groups and hierarchy ----------------

# Author: Kate P Maia
# Checked: 12/2023

# Reads data generated in previous scripts to plot manuscript figures on group structure of ecological networks and hierarchical structure. 

# ESTOU USANDO O N certo em todos os plots?

# --------------------- Loading library, code and data --------------------

library(tidyverse)
library(RColorBrewer)
library(reshape2)
library(ghibli)

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

# multLC & starLC (to be removed plot by plot)
multLC <- c("Herb_CuevasReyes", "IZ_PA_DIM100", "IZ_PA_PAcamp", "M_AF_002_02", "M_AF_002_11")
starLC <- c("IZ_PA_DIM10", "IZ_PA_PA10", "M_AF_002_07", "M_AF_002_16", "M_PL_061_33")

dim(net_struct); length(unique(sp_struct$ID))

pal <- c("#D85C00", "#31489F")

# -------------------------------------------------------------------------
# ----------- Fig 2: Plotting module-sector congruence results ------------

data <- net_struct %>% filter(!ID %in% c(multLC, starLC)) %>% 
  select(ID, IntSign, Code, MS_congr, MS_mean, MS_sd, MS_pval)

pa <- data %>% ggplot(aes(x = Code, y = MS_congr, fill = IntSign)) +
  geom_boxplot(color = c(rep(pal[1], 4), rep(pal[2], 4)), alpha = 0.5) +
  xlab("") + ylab("Average module-sector congruence") +
  coord_cartesian(ylim = c(0.50, 1.00)) + scale_fill_manual(values = pal) + 
  theme_classic() + theme(axis.text.x = element_text(size = 14),axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 14), legend.position = "none"); pa

pb <- data %>% select(Code, MS_pval) %>% 
  mutate(plogic = ifelse(MS_pval <= 0.05, TRUE, FALSE)) %>% 
  mutate(plogic = factor(plogic, levels = c(FALSE, TRUE))) %>% 
  group_by(Code, plogic) %>% summarise(count = n()) %>% 
  ggplot(aes(x = Code, y = count, fill = Code, alpha = plogic)) + 
  geom_bar(position = "fill", stat = "identity") + # position = "stack"
  xlab("") + ylab("Proportion of networks") + scale_alpha_manual(values = c(0.2, 0.9)) + 
  scale_fill_manual(values = c(rep(pal[1], 4), rep(pal[2], 4))) +
  ggtitle("Module-sector congruence: null model") +
  theme_classic() + theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14), plot.title = element_text(size = 16), legend.position = "none")

cowplot::plot_grid(pa, pb, labels = c("a)", "b)"), nrow = 1, ncol = 2, label_size = 12, label_fontfamily = "", label_fontface = "plain")

# -------------------------------------------------------------------------
# ------------- Fig 3: Plotting subgraphs in modules/sectors -------------- 

pa <- net_struct %>% filter(!ID %in% c(multLC, starLC)) %>% 
  select(ID, IntSign, Code, n_m6, n_m6_sect) %>%
  mutate(In = n_m6_sect/n_m6) %>% filter(!is.nan(In)) %>% # removes 21 with no m6
  ggplot(aes(x = Code, y = In, col = IntSign, fill = IntSign)) + geom_boxplot(alpha = 0.5) +
  scale_color_manual(values = pal) + scale_fill_manual(values = pal) + 
  xlab("") + ylab("Proportion of subgraphs inside sectors") +
  theme_classic() + theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14), plot.title = element_text(size = 16), legend.position = "none"); pa

pb <- net_struct %>% filter(!ID %in% c(multLC, starLC)) %>% 
  select(ID, IntSign, Code, p_m6_sect) %>%
  mutate(p_m6_sect = ifelse(p_m6_sect <= 0.05, TRUE, FALSE)) %>% 
  mutate(p_m6_sect = factor(p_m6_sect, levels = c(FALSE, TRUE))) %>% 
  group_by(Code, p_m6_sect) %>% summarise(count = n()) %>% 
  ggplot(aes(x = Code, y = count, fill = Code, alpha = p_m6_sect)) + 
  geom_bar(position = "fill", stat = "identity") + # position = "stack"
  xlab("") + ylab("Proportion of networks") + ggtitle("Subgraphs inside sectors: null model") +
  scale_alpha_manual(values = c(0.2, 0.9)) + scale_fill_manual(values = rep(pal, each = 4)) +
  theme_classic() + theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14), plot.title = element_text(size = 16), legend.position = "none"); pb

pc <- net_struct %>%  filter(!ID %in% c(multLC, starLC)) %>% 
  select(ID, IntSign, Code, n_m6, n_m6_mod) %>%
  mutate(In = n_m6_mod/n_m6) %>% filter(!is.nan(In)) %>% # NaNs are 0/0
  ggplot(aes(x = Code, y = In, col = IntSign, fill = IntSign)) + geom_boxplot(alpha = 0.5) +
  scale_color_manual(values = pal) + scale_fill_manual(values = pal) + 
  xlab("") + ylab("Proportion of subgraphs inside modules") + theme_classic() + 
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14), plot.title = element_text(size = 16), legend.position = "none"); pc

pd <- net_struct %>%  filter(!ID %in% c(multLC, starLC)) %>% 
  select(ID, IntSign, Code, p_m6_mod) %>%
  mutate(p_m6_mod = ifelse(p_m6_mod <= 0.05, TRUE, FALSE)) %>% 
  mutate(p_m6_mod = factor(p_m6_mod, levels = c(FALSE, TRUE))) %>% 
  group_by(Code, p_m6_mod) %>% summarise(count = n()) %>% 
  ggplot(aes(x = Code, y = count, fill = Code, alpha = p_m6_mod)) + 
  geom_bar(position = "fill", stat = "identity") + # position = "stack"
  xlab("") + ylab("Proportion of networks") + ggtitle("Subgraphs inside modules: null model") +
  scale_alpha_manual(values = c(0.2, 0.9)) + scale_fill_manual(values = rep(pal, each = 4)) +
  theme_classic() + theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14), plot.title = element_text(size = 16), legend.position = "none"); pd

cowplot::plot_grid(pa, pb, pc, pd, labels = c("a)", "b)", "c)", "d)"), nrow = 2, ncol = 2, label_size = 14, label_fontfamily = "", label_fontface = "plain")

# -------------------------------------------------------------------------
# ------------- Fig S1: Conceptual figure (use of no data)  ---------------

# -------------------------------------------------------------------------
# ----------- Fig S2: Plotting component size x no components -------------

net_struct %>% select(ID, IntSign, IntType, 2:8) %>% 
  mutate(PropLarge = LargeCompSize/(NRow + NCol)) %>% 
  ggplot(aes(x = NComp, y = PropLarge, colour = IntSign)) + geom_point(alpha = 0.5) +
  xlab("Number of network components") + ylab("Proportion of species in the largest component") + theme_classic() + 
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 14), axis.title.x = element_text(size = 14), strip.text = element_text(size = 12), legend.position = "none") + 
  scale_colour_manual(values = pal) + facet_wrap(~IntType, nrow = 2)

# -------------------------------------------------------------------------
# --------------- Fig S3: Plotting balance in sector size -----------------

ponyo <- ghibli_palette("PonyoMedium")[1:7]
harry <- harrypotter::hp(18, house = "LunaLovegood")

net_struct <- net_struct %>% filter(!ID %in% c(multLC, starLC)) 
sp_struct <- sp_struct %>% filter(!ID %in% c(multLC, starLC)) 

sector_tab <- sp_struct %>% filter(!is.na(GC_Aff)) %>% # # 28930 sp in GC1 
  group_by(ID) %>% count(Fiedler) %>% 
  pivot_wider(names_from = Fiedler, values_from = n) %>% 
  mutate(Sct1 = `-1` / (`-1` + `1`), Sct2 = `1` / (`-1` + `1`))
all(rowSums(sector_tab[,2:3]) == net_struct[!net_struct$ID %in% c(multLC, starLC), "LargeCompSize"])
all(rowSums(sector_tab[,4:5]) == 1)

sector_tab$Code <- net_struct[match(sector_tab$ID, net_struct$ID), "Code"]

sector_tab <- sector_tab %>% select(-c(2:3)) %>% 
  pivot_longer(2:3, names_to = "Sector", values_to = "Size")

pa <- sector_tab %>% group_by(Code, Sector) %>% summarise(Size = mean(Size)) %>% 
  ggplot(aes(x = Code, y = Size, fill = Sector)) +
  geom_bar(position = "fill", stat = "identity") + 
  labs(x = "", y = "Proportion of species") + 
  scale_fill_manual(values = c(ponyo[5], ponyo[6])) + theme_classic() +
  theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 14)); pa

data <- sector_tab %>% filter(Sector == "Sct1") %>% 
  mutate(Size = ifelse(Size > 0.5, 1 - Size, Size)) # so all Sizes are below 0.5
data$Balance <- ifelse(data$Size <= 0.15, "Low", data$Size)
data$Balance <- ifelse(data$Size > 0.15 & data$Size <= 0.35, "Medium", data$Balance)
data$Balance <- ifelse(data$Size > 0.35 & data$Size <= 0.5, "High", data$Balance)
data <- data %>% mutate(Balance = factor(Balance, levels = c( "High", "Medium", "Low")))
data %>% group_by(Balance) %>% summarise(min = min(Size), max = max(Size))
data %>% group_by(Balance) %>% count() # 366 = 128 + 152 + 86

pb <- data %>% group_by(Code, Balance) %>% count(Code) %>% 
  ggplot(aes(x = Code, y = n, fill = Balance)) +
  geom_bar(position = "fill", stat = "identity") + 
  labs(x = "", y = "Proportion of networks") +
  scale_fill_manual(values = harry[c(3,6,9)]) + theme_classic() +
  theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 14)); pb

cowplot::plot_grid(pa, pb, labels = c("a)", "b)"), nrow = 1, label_size = 12, label_fontfamily = "", label_fontface = "plain")

# -------------------------------------------------------------------------
# ------- Fig SX: Plotting correlation between modularity indexes ---------

data <- net_struct %>% select(ID, SA_NMod, SA_Modularity, FG_NMod, FG_Modularity)

pa <- ggplot(data, aes(x = SA_Modularity, y = FG_Modularity)) + geom_jitter(alpha = 0.3) + 
  xlab(expression(paste("Simulated Annealing - ", Q[B]))) + ylab("Fast Greedy - Q") + 
  ggtitle("Modularity") + theme_classic() +
  theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14), plot.title = element_text(size = 16))

pb <- ggplot(data, aes(x = SA_NMod, y = FG_NMod)) + geom_jitter(alpha = 0.3) + 
  xlab(expression(paste("Simulated Annealing - ", Q[B]))) + ylab("Fast Greedy - Q") + 
  ggtitle("Number of modules") + theme_classic() +
  theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14), plot.title = element_text(size = 16))

cowplot::plot_grid(pa, pb, labels = c("a)", "b)"), nrow = 1, label_size = 12, label_fontfamily = "", label_fontface = "plain")

# -------------------------------------------------------------------------
# ------------- Fig S4: Plotting network modularity results ---------------

data <- net_struct %>% select(ID, IntSign, Code, Connectance, SA_Modularity, FG_PNull)

pa <- data %>% ggplot(aes(x = IntSign, y = SA_Modularity)) + 
  geom_violin() + geom_boxplot(width = 0.1, color = pal, fill = pal, alpha = 0.2) +
  xlab("") + theme_classic() + 
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 14)); pa

pb <- data %>% ggplot(aes(x = Connectance, y = SA_Modularity, col = IntSign)) +
  geom_point(alpha = 0.5) + scale_color_manual(values = pal) + theme_classic() + 
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 14), axis.title.x = element_text(size = 14), legend.position = "none"); pb

pc <- data %>% ggplot(aes(x = Code, y = SA_Modularity)) + 
  geom_boxplot(colour = c(rep(pal[1], 4), rep(pal[2], 4)), fill = c(rep(pal[1], 4), rep(pal[2], 4)), alpha = 0.5) +
  xlab("") + theme_classic() +
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 14)); pc

pd <- data %>% mutate(FG_PNull = ifelse(FG_PNull <= 0.05, TRUE, FALSE)) %>% 
  select(ID, FG_PNull, Code) %>% group_by(Code) %>% count(FG_PNull) %>% 
  ggplot(aes(x = Code, y = n, fill = Code, alpha = FG_PNull)) +
  geom_bar(position = "fill", stat = "identity") + xlab("") +
  ylab("Proportion of modular networks") + 
  scale_alpha_manual(values = c(0.2, 0.9)) +
  scale_fill_manual(values = rep(pal, each = 4)) +
  theme_classic() + theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 14), axis.title.x = element_text(size = 14), legend.position = "none"); pd

cowplot::plot_grid(pa, pb, pc, pd, labels = c("a)", "b)", "c)", "d)"), nrow = 2, label_size = 12, label_fontfamily = "", label_fontface = "plain")

# -------------------------------------------------------------------------
# -------------- Fig S5: Plotting no subgraphs x net size -----------------

net_struct %>% filter(n_m6 != 0) %>% # removes 21 with no m6
  ggplot(aes(x = log(NRow + NCol), y = log(n_m6), color = IntType)) + geom_point(alpha = 0.5) +
  xlab("Network size (log)") + ylab("Number of subgraphs (log)") +
  scale_color_manual(values = rep(pal, each = 4)) + theme_classic() +
  theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14), strip.text = element_text(size = 12), legend.position = "none") + facet_wrap(~IntType, nrow = 2)

# -------------------------------------------------------------------------
# -------------- Fig S6: Plotting no subgraphs x net size -----------------

pa <- net_struct %>% select(Code, ID, n_m6, n_m6_ms) %>% # all network types
  filter(!is.na(n_m6_ms)) %>% filter(n_m6 != 0) %>% # removes M_PL_062, 21 no m6 (see above)
  mutate(Sign = ifelse(Code %in% c("Pa-Ho", "Ph-B", "Her-P", "Pr-Pr"), "A", "M")) %>%
  mutate(In = n_m6_ms/n_m6) %>% filter(!is.na(In)) %>% # not NA left to remove 
  ggplot(aes(x = Code, y = In, col = Sign, fill = Sign)) + geom_boxplot(alpha = 0.5) +
  scale_color_manual(values = pal) + scale_fill_manual(values = pal) +
  xlab("") + ylab("Proportion of subgraphs inside sectors and modules") + theme_classic() +
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 12), legend.position = "none"); pa

pb <- net_struct %>% filter(!is.na(p_m6_ms)) %>%  # removes M_PL_062 (see above) 
  mutate(p_m6_ms = ifelse(p_m6_ms <= 0.05, TRUE, FALSE)) %>% 
  mutate(p_m6_ms = factor(p_m6_ms, levels = c(FALSE, TRUE))) %>% 
  group_by(Code, p_m6_ms) %>% summarise(count = n()) %>% 
  ggplot(aes(x = Code, y = count, fill = Code, alpha = p_m6_ms)) + 
  geom_bar(position = "fill", stat = "identity") + # position = "stack"
  xlab("") + ylab("Proportion of networks") +
  scale_alpha_manual(values = c(0.2, 0.9)) + scale_fill_manual(values = rep(pal, each = 4)) +
  ggtitle("Subgraphs inside sectors and modules: null model") + theme_classic() +
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 12), plot.title = element_text(size = 14), legend.position = "none"); pb

cowplot::plot_grid(pa, pb, labels = c("a)", "b)"), ncol = 2, label_size = 14, label_fontfamily = "", label_fontface = "plain")

# -------------------------------------------------------------------------
# ------------- Fig S7: Plotting proportion of in-group links -------------

head(net_struct)

net_struct %>% select(ID, IntSign, IntType, LinkSct, LinkMod, LinkMot) %>% 
  rename(Sct = LinkSct, Mod = LinkMod, Subg = LinkMot) %>% 
  pivot_longer(4:6, names_to = "Group", values_to = "Value") %>% 
  mutate(Group = factor(Group, levels = c("Subg", "Mod", "Sct"))) %>%
  ggplot(aes(x = Group, y = Value, color = IntSign, fill = IntSign)) + geom_boxplot(alpha = 0.5) +
  xlab("") + ylab("Proportion of links inside groups") +
  facet_wrap(~IntType, nrow = 2) + theme_classic() +
  scale_color_manual(values = pal) + scale_fill_manual(values = pal) + 
  theme(legend.position = "none", axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 14), strip.text = element_text(size = 12))