# --------------- Manuscript figures: groups and hierarchy ----------------

# Author: Kate P Maia
# Checked: 08/2024

# Reads data generated in previous scripts to plot manuscript figures on group structure of ecological networks and hierarchical structure. 

# --------------------- Loading library, code and data --------------------

library(tidyverse)
library(RColorBrewer)
library(cowplot)
library(reshape2)
library(ghibli)
source("./Scripts/Functions/zscore.R")

net_struct <- read.table("./Outputs/05_Net-Level_Struct.txt", header = T, stringsAsFactors = F)
sp_struct <- read.table("./Outputs/04_Node-Level_Struct.txt", header = T, stringsAsFactors = F)

sp_struct <- sp_struct[order(sp_struct$ID),]; all(net_struct$ID == levels(sp_struct$ID))

net_struct$IntType <- factor(net_struct$IntType, levels = c("Host-Parasite", "Phage-Bacteria", "Plant-Herbivore", "Food Webs", "Pollination", "Seed Dispersal", "Anemone-Fish", "Plant-Ant"))
net_struct$Code <- net_struct$IntSign <- net_struct$IntType
levels(net_struct$IntType) <- c("Parasite-Host", "Phage-Bacteria", "Herbivore-Plant", "Predator-Prey", "Pollinator-Plant", "Seed Disperser-Plant", "Fish-Anemone", "Ant-Plant")
levels(net_struct$IntSign) <- rep(c("A", "M"), each = 4)
levels(net_struct$Code) <- c("Pa-Ho", "Ph-B", "Her-P", "Pr-Pr", "Pol-P", "SD-P", "F-Ane", "Ant-P")

sp_struct$IntType <- net_struct[match(sp_struct$ID, net_struct$ID), "IntType"]
sp_struct$IntSign <- net_struct[match(sp_struct$ID, net_struct$ID), "IntSign"]

# multLC & starLC (to be removed plot by plot)
multLC <- c("Herb_CuevasReyes", "IZ_PA_DIM100", "IZ_PA_PAcamp", "M_AF_002_02", "M_AF_002_11")
starLC <- c("IZ_PA_DIM10", "IZ_PA_PA10", "M_AF_002_07", "M_AF_002_16", "M_PL_061_33")

pal <- c("#D85C00", "#31489F")

# -------------------------------------------------------------------------
# --------------- Congruence: zscore and normalised zscore ----------------

temp <- net_struct %>% filter(!ID %in% c(multLC, starLC)) # 366 ID to calculate zscore
temp <- temp %>% select(ID, IntSign, Code, 15:31); head(temp); dim(temp) # facilitate visualization
temp %>% filter(n_m6 == 0) # 21 networks with no m6

# module-sector congruence
msc <- temp %>% select(ID, Code, MS_congr, MS_mean, MS_sd); head(msc)
msc <- zscore(msc); msc %>% filter(is.na(Z) | is.na(NormZ))

temp <- temp %>% filter(n_m6 != 0); dim(temp) # if no subgraphs, we do not test for hierarchy 

# subgraph-sector congruence
m6sc <- temp %>% select(ID, Code, n_m6_sect, mean_m6_sect, sd_m6_sect); head(m6sc)
m6sc %>% filter(sd_m6_sect == 0) # 6 0 in sd congruence
m6sc <- zscore(m6sc); m6sc %>% filter(is.na(Z)) # 6 NaN (0/0)

# subgraph-module congruence
m6mc <- temp %>% select(ID, Code, n_m6_mod, mean_m6_mod, sd_m6_mod); head(m6mc)
m6mc %>% filter(sd_m6_mod == 0) # 24 0 in sd congruence
m6mc <- zscore(m6mc); m6mc %>% filter(is.na(Z)) # 24 Sd = 0

# subgraph-module-sector congruence ("double hierarchy")
m6msc <- temp %>% filter(ID != "M_PL_062") %>% # additionally removes M_PL_062 - too large
  select(ID, Code, n_m6_ms, mean_m6_ms, sd_m6_ms); head(m6msc); dim(m6msc) # 344
m6msc %>% filter(sd_m6_ms == 0) # 38 0 in sd congruence
m6msc <- zscore(m6msc); m6msc %>% filter(is.na(Z)) # 36 Sd = 0
dim(m6msc) # 2 removed due to Inf

zdata <- left_join(msc, m6sc, by = c("ID", "Code"), suffix = c("_ms", "_m6s")) %>% 
  left_join(m6mc, by = c("ID", "Code")) %>% 
  left_join(m6msc, by = c("ID", "Code"), suffix = c("_m6m", "_m6ms")); dim(zdata)

# -------------------------------------------------------------------------
# ----------- Fig 2: Plotting module-sector congruence results ------------

pa <- net_struct %>% filter(!ID %in% c(multLC, starLC)) %>% # N = 366 
  select(ID, IntSign, Code, contains("MS_")) %>% 
  ggplot(aes(x = Code, y = MS_congr, fill = IntSign)) +
  geom_boxplot(color = c(rep(pal[1], 4), rep(pal[2], 4)), alpha = 0.5) +
  xlab("") + ylab("Module-sector congruence") +
  coord_cartesian(ylim = c(0.00, 1.00)) + scale_fill_manual(values = pal) + 
  theme_classic() + theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 14), legend.position = "none"); pa

label <- net_struct %>% filter(!ID %in% c(multLC, starLC)) %>% group_by(Code) %>% count()
pb <- net_struct %>% filter(!ID %in% c(multLC, starLC)) %>% # N = 366
  select(ID, IntSign, Code, MS_pval) %>% # N = 366 
  mutate(plogic = ifelse(MS_pval <= 0.05, TRUE, FALSE)) %>% 
  mutate(plogic = factor(plogic, levels = c(FALSE, TRUE))) %>% 
  group_by(Code, plogic) %>% summarise(count = n()) %>% 
  ggplot(aes(x = Code, y = count, fill = Code, alpha = plogic)) + 
  geom_bar(position = "fill", stat = "identity") + # position = "stack"
  xlab("") + ylab("Proportion of networks") + scale_alpha_manual(values = c(0.2, 0.9)) + 
  scale_fill_manual(values = c(rep(pal[1], 4), rep(pal[2], 4))) +
  ggtitle("Module-sector congruence: null model") +
  theme_classic() + theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 14), plot.title = element_text(size = 14), legend.position = "none")
pb <- pb + annotate("text", x = 1:8, y = 1.02, label = label$n); pb

pc <- zdata %>% ggplot() + geom_histogram(aes(x = Z_ms), alpha = 0.6) + # 366
  xlab("Module-sector congruence (z-score)") + ylab("Frequecy") + theme_cowplot(); pc

cowplot::plot_grid(pa, pb, pc, labels = c("a)", "b)", "c)"), nrow = 1, label_size = 12, label_fontfamily = "", label_fontface = "plain")

# -------------------------------------------------------------------------
# ------------- Fig 3: Plotting subgraphs in modules/sectors -------------- 

pa <- net_struct %>% filter(!ID %in% c(multLC, starLC)) %>% # N = 366
  filter(n_m6 != 0) %>% # N = 345 (366 - 21 with 0 m6)
  select(ID, IntSign, Code, n_m6, n_m6_sect) %>%
  mutate(In = n_m6_sect/n_m6) %>%
  ggplot(aes(x = Code, y = In, col = IntSign, fill = IntSign)) + geom_boxplot(alpha = 0.5) +
  scale_color_manual(values = pal) + scale_fill_manual(values = pal) + 
  xlab("") + ylab("Subgraph-sector congruence") +
  theme_classic() + theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 14), legend.position = "none"); pa

label <- net_struct %>% filter(!ID %in% c(multLC, starLC)) %>% filter(n_m6 != 0) %>% group_by(Code) %>% count()
pb <- net_struct %>% filter(!ID %in% c(multLC, starLC)) %>% # N = 366
  filter(n_m6 != 0) %>% # N = 345 (366 - 21 with 0 m6)
  mutate(p_m6_sect = ifelse(p_m6_sect <= 0.05, TRUE, FALSE)) %>% 
  mutate(p_m6_sect = factor(p_m6_sect, levels = c(FALSE, TRUE))) %>% 
  group_by(Code, p_m6_sect) %>% summarise(count = n()) %>% # 15/16 potential groups
  ggplot(aes(x = Code, y = count, fill = Code, alpha = p_m6_sect)) + 
  geom_bar(position = "fill", stat = "identity") + # position = "stack"
  xlab("") + ylab("Proportion of networks") + ggtitle("Subgraph-sector congruence: null model") +
  scale_alpha_manual(values = c(0.2, 0.9)) + scale_fill_manual(values = rep(pal, each = 4)) +
  theme_classic() + theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 14), plot.title = element_text(size = 14), legend.position = "none")
pb <- pb + annotate("text", x = 1:8, y = 1.02, label = label$n); pb

# 27 removed: 21 no m6 (366 - 345) + 6 sd = 0 (see above)
pc <- zdata %>% ggplot() + geom_histogram(aes(x = Z_m6s), alpha = 0.6) + # N = 339
  xlab("Subgraph-sector congruence (z-score)") + ylab("Frequecy") + 
  theme_cowplot() + theme(axis.title.x = element_text(size = 12)); pc

pd <- net_struct %>% filter(!ID %in% c(multLC, starLC)) %>% # N = 366
  filter(n_m6 != 0) %>% # N = 345 (366 - 21 with 0 m6)
  select(ID, IntSign, Code, n_m6, n_m6_mod) %>% # N = 345 (366 - 21 with 0 m6)
  mutate(In = n_m6_mod/n_m6) %>%
  ggplot(aes(x = Code, y = In, col = IntSign, fill = IntSign)) + geom_boxplot(alpha = 0.5) +
  scale_color_manual(values = pal) + scale_fill_manual(values = pal) + 
  xlab("") + ylab("Subgraph-module congruence") + theme_classic() + 
  theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 14), legend.position = "none"); pd

pe <- net_struct %>% filter(!ID %in% c(multLC, starLC)) %>% # N = 366
  filter(n_m6 != 0) %>% # N = 345 (366 - 21 with 0 m6)
  select(ID, IntSign, Code, p_m6_mod) %>%
  mutate(p_m6_mod = ifelse(p_m6_mod <= 0.05, TRUE, FALSE)) %>% 
  mutate(p_m6_mod = factor(p_m6_mod, levels = c(FALSE, TRUE))) %>% 
  group_by(Code, p_m6_mod) %>% summarise(count = n()) %>% 
  ggplot(aes(x = Code, y = count, fill = Code, alpha = p_m6_mod)) + 
  geom_bar(position = "fill", stat = "identity") + # position = "stack"
  xlab("") + ylab("Proportion of networks") + ggtitle("Subgraph-module congruence: null model") +
  scale_alpha_manual(values = c(0.2, 0.9)) + scale_fill_manual(values = rep(pal, each = 4)) +
  theme_classic() + theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 14), plot.title = element_text(size = 14), legend.position = "none"); pe

# 45 removed: 21 no m6 (366 - 345) + 24 sd = 0 (see above)
pf <- zdata %>% ggplot() + geom_histogram(aes(x = Z_m6m), alpha = 0.6) + 
  xlab("Subgraph-module congruence (z-score)") + ylab("Frequecy") + 
  theme_cowplot() + theme(axis.title.x = element_text(size = 12)); pf

cowplot::plot_grid(pa, pb, pc, pd, pe, pf, labels = c("a)", "b)", "c)", "d)", "e)", "f)"), nrow = 2, label_size = 14, label_fontfamily = "", label_fontface = "plain")

# -------------------------------------------------------------------------
# ---------------- Plotting component size x no components ----------------

net_struct %>% select(ID, IntSign, IntType, 2:8) %>% 
  mutate(PropLarge = LargeCompSize/(NRow + NCol)) %>% 
  ggplot(aes(x = NComp, y = PropLarge, colour = IntSign)) + geom_point(alpha = 0.5) +
  xlab("Number of network components") + ylab("Proportion of species in the largest component") + theme_classic() + 
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 14), axis.title.x = element_text(size = 14), strip.text = element_text(size = 12), legend.position = "none") + 
  scale_colour_manual(values = pal) + facet_wrap(~IntType, nrow = 2)

# -------------------------------------------------------------------------
# -------------------- Plotting balance in sector size --------------------

ponyo <- ghibli_palette("PonyoMedium")[1:7]
harry <- harrypotter::hp(18, house = "LunaLovegood")

sector_tab <- sp_struct %>% filter(!ID %in% c(multLC, starLC)) %>% 
  filter(!is.na(GC_Aff)) %>% # # 28930 sp in GC1 
  group_by(ID) %>% count(Fiedler) %>% pivot_wider(names_from = Fiedler, values_from = n) %>% 
  mutate(Sct1 = `-1` / (`-1` + `1`), Sct2 = `1` / (`-1` + `1`)) %>% select(-c(2:3)) %>% 
  mutate(Large = ifelse(Sct1 >= Sct2, Sct1, Sct2), Small = ifelse(Sct1 >= Sct2, Sct2, Sct1))
sector_tab$Code <- net_struct[match(sector_tab$ID, net_struct$ID), "Code"]

pa <- sector_tab %>% pivot_longer(4:5, names_to = "Sector", values_to = "Size") %>% 
  group_by(Code, Sector) %>% summarise(Size = mean(Size)) %>% 
  ggplot(aes(x = Code, y = Size, fill = Sector)) +
  geom_bar(position = "fill", stat = "identity") + 
  labs(x = "", y = "Proportion of species") + 
  scale_fill_manual(values = c(ponyo[5], ponyo[6])) + theme_classic() +
  theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 14)); pa

sector_tab$Balance <- ifelse(sector_tab$Small <= 0.15, "Low", sector_tab$Small)
sector_tab$Balance <- ifelse(sector_tab$Small > 0.15 & sector_tab$Small <= 0.35, "Medium", sector_tab$Balance)
sector_tab$Balance <- ifelse(sector_tab$Small > 0.35 & sector_tab$Small <= 0.5, "High", sector_tab$Balance)
sector_tab <- sector_tab %>% mutate(Balance = factor(Balance, levels = c("High", "Medium", "Low")))
sector_tab %>% group_by(Balance) %>% summarise(min = min(Small), max = max(Small))

pb <- sector_tab %>% group_by(Code, Balance) %>% count(Code) %>% 
  ggplot(aes(x = Code, y = n, fill = Balance)) +
  geom_bar(position = "fill", stat = "identity") + 
  labs(x = "", y = "Proportion of networks") +
  scale_fill_manual(values = harry[c(3,6,9)]) + theme_classic() +
  theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 14)); pb

cowplot::plot_grid(pa, pb, labels = c("a)", "b)"), nrow = 1, label_size = 12, label_fontfamily = "", label_fontface = "plain")

# -------------------------------------------------------------------------
# ------------ Plotting correlation between modularity indexes ------------

data <- net_struct %>% filter(!ID %in% c(multLC, starLC)) %>% # N = 366
  select(ID, SA_NMod, SA_Modularity, FG_NMod, FG_Modularity)

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
# ------------------ Plotting network modularity results ------------------

data <- net_struct %>% filter(!ID %in% c(multLC, starLC)) %>% # N = 366
  select(ID, IntSign, Code, Connectance, SA_Modularity, FG_PNull)

pa <- data %>% mutate(IntSign = recode_factor(IntSign, A = "Antagonistic", M = "Mutualistic")) %>% 
  ggplot(aes(x = IntSign, y = SA_Modularity)) + 
  geom_violin() + geom_boxplot(width = 0.1, color = pal, fill = pal, alpha = 0.2) +
  xlab("") + ylab("Modularity") + theme_classic() + 
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 14)); pa

pb <- data %>% ggplot(aes(x = Connectance, y = SA_Modularity, col = IntSign)) +
  geom_point(alpha = 0.5) + scale_color_manual(values = pal) + 
  ylab("Modularity") + theme_classic() + 
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 14), axis.title.x = element_text(size = 14), legend.position = "none"); pb

pc <- data %>% ggplot(aes(x = Code, y = SA_Modularity)) + 
  geom_boxplot(colour = c(rep(pal[1], 4), rep(pal[2], 4)), fill = c(rep(pal[1], 4), rep(pal[2], 4)), alpha = 0.5) +
  xlab("") + ylab("Modularity") + theme_classic() +
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
# ------------------- Plotting no subgraphs x net size --------------------

net_struct %>% filter(!ID %in% c(multLC, starLC)) %>% # N = 366
  filter(n_m6 != 0) %>% # removes 21 with no m6
  ggplot(aes(x = log(NRow + NCol), y = log(n_m6), color = IntType)) + geom_point(alpha = 0.5) +
  xlab("Network size (log)") + ylab("Number of subgraphs (log)") +
  scale_color_manual(values = rep(pal, each = 4)) + theme_classic() +
  theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14), strip.text = element_text(size = 12), legend.position = "none") + facet_wrap(~IntType, nrow = 2)

# -------------------------------------------------------------------------
# ------------------- Subgraph-Module-Sector Congruence -------------------

pa <- net_struct %>% filter(!ID %in% c(multLC, starLC)) %>% # N = 366
  filter(n_m6 != 0) %>% filter(ID != "M_PL_062") %>% # N = 344 (366 - 21 with 0 m6)
  select(ID, Code, IntSign, n_m6, n_m6_ms) %>% mutate(In = n_m6_ms/n_m6) %>%
  ggplot(aes(x = Code, y = In, col = IntSign, fill = IntSign)) + geom_boxplot(alpha = 0.5) +
  scale_color_manual(values = pal) + scale_fill_manual(values = pal) +
  xlab("") + ylab("Subgraphs-module-sector congruence") + theme_classic() +
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 12), legend.position = "none"); pa

label <- net_struct %>% filter(!ID %in% c(multLC, starLC)) %>% filter(n_m6 != 0) %>% filter(ID != "M_PL_062") %>% group_by(Code) %>% count()
pb <- net_struct %>% filter(!ID %in% c(multLC, starLC)) %>% # N = 366
  filter(n_m6 != 0) %>% filter(ID != "M_PL_062") %>% # N = 344 (366 - 21 with 0 m6)
  mutate(p_m6_ms = ifelse(p_m6_ms <= 0.05, TRUE, FALSE)) %>% 
  mutate(p_m6_ms = factor(p_m6_ms, levels = c(FALSE, TRUE))) %>% 
  group_by(Code, p_m6_ms) %>% summarise(count = n()) %>% 
  ggplot(aes(x = Code, y = count, fill = Code, alpha = p_m6_ms)) + 
  geom_bar(position = "fill", stat = "identity") + # position = "stack"
  xlab("") + ylab("Proportion of networks") +
  scale_alpha_manual(values = c(0.2, 0.9)) + scale_fill_manual(values = rep(pal, each = 4)) +
  ggtitle("Subgraphs-module-sector congruence: null model") + theme_classic() +
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 12), plot.title = element_text(size = 14), legend.position = "none")
pb <- pb + annotate("text", x = 1:8, y = 1.02, label = label$n); pb

cowplot::plot_grid(pa, pb, labels = c("a)", "b)"), ncol = 2, label_size = 14, label_fontfamily = "", label_fontface = "plain")

# -------------------------------------------------------------------------
# ----------------- Plotting proportion of in-group links -----------------

net_struct %>% filter(!ID %in% c(multLC, starLC)) %>%
  select(ID, IntSign, IntType, LinkSct, LinkMod, LinkMot) %>% 
  rename(Sct = LinkSct, Mod = LinkMod, Subg = LinkMot) %>% 
  pivot_longer(4:6, names_to = "Group", values_to = "Value") %>% # NA for M_PL_062
  mutate(Group = factor(Group, levels = c("Subg", "Mod", "Sct"))) %>%
  ggplot(aes(x = Group, y = Value, color = IntSign, fill = IntSign)) + geom_boxplot(alpha = 0.5) +
  xlab("") + ylab("Proportion of links inside groups") +
  facet_wrap(~IntType, nrow = 2) + theme_classic() +
  scale_color_manual(values = pal) + scale_fill_manual(values = pal) + 
  theme(legend.position = "none", axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 14), strip.text = element_text(size = 12))

# -------------------------------------------------------------------------
# --------------- Normalised zscore congruence histograms  ----------------

pa <- zdata %>% ggplot() + geom_histogram(aes(x = NormZ_ms), alpha = 0.6) + 
  xlab("Module-sector congruence \n (normalised z-score)") + 
  ylab("Frequecy") + theme_cowplot(); pa

pb <- zdata %>% ggplot() + geom_histogram(aes(x = NormZ_m6s), alpha = 0.6) + 
  xlab("Subgraph-sector congruence \n (normalised z-score)") + 
  ylab("Frequecy") + theme_cowplot(); pb

pc <- zdata %>% ggplot() + geom_histogram(aes(x = NormZ_m6m), alpha = 0.6) + 
  xlab("Subgraph-module congruence \n (normalised z-score)") + 
  ylab("Frequecy") + theme_cowplot(); pc

cowplot::plot_grid(pa, pb, pc, labels = c("a)", "b)", "c)"), nrow = 1, label_size = 14, label_fontfamily = "", label_fontface = "plain")

# -------------------------------------------------------------------------
# -------------------- Congruence Significance Profile --------------------

# data to explore missing values
miss <- zdata %>% mutate(Z_ms = is.na(Z_ms), NormZ_ms = is.na(NormZ_ms), Z_m6s = is.na(Z_m6s), NormZ_m6s = is.na(NormZ_m6s), Z_m6m = is.na(Z_m6m), NormZ_m6m = is.na(NormZ_m6m), Z_m6ms = is.na(Z_m6ms), NormZ_m6ms = is.na(NormZ_m6ms)) %>% 
  pivot_longer(3:10, names_to = "zscore", values_to = "value") %>% 
  group_by(Code, zscore, value) %>% count()

# ----- ZSCORE ----- #

A <- zdata %>% filter(Code %in% c("Pa-Ho", "Ph-B", "Her-P", "Pr-Pr")) %>%
  select(ID, Code, !contains("Norm")) %>%
  pivot_longer(3:6, names_to = "cngr", values_to = "value") %>% 
  mutate(cngr = factor(cngr, levels = c("Z_ms", "Z_m6s", "Z_m6m", "Z_m6ms"))) %>% 
  mutate(cngr = recode_factor(cngr, "Z_ms" = "1", "Z_m6s" = "2", "Z_m6m" = "3", "Z_m6ms" = "4")) %>% # trick: numeric to plot geom_line
  mutate(cngr = as.numeric(cngr)) %>% 
  ggplot(aes(x = cngr, y = value, color = ID)) + 
  geom_line(alpha = 0.2) + geom_point(alpha = 0.6) + 
  scale_color_manual(values = rep(pal[1], 154)) + 
  scale_x_continuous(labels = c("Mod-Sct", "Subg-Sct", "Subg-Mod", "Subg-\nMod-Sct")) + 
  xlab("") + ylab("Congruence coefficients (z-scores)") +
  facet_wrap(~Code, nrow = 4, scales = "free_y") + theme_cowplot() + 
  theme(legend.position = "none", plot.margin = unit(c(0,0,0,0), "cm")); A

miss %>% filter(Code %in% c("Pa-Ho", "Ph-B", "Her-P", "Pr-Pr")) %>% 
  filter(zscore %in% c("Z_ms", "Z_m6s", "Z_m6m", "Z_m6ms")) %>% 
  filter(value == T) %>% ungroup() %>% select(n) %>% sum()

M <- zdata %>% filter(Code %in% c("Pol-P", "SD-P", "F-Ane", "Ant-P")) %>%
  select(ID, Code, !contains("Norm")) %>%
  pivot_longer(3:6, names_to = "cngr", values_to = "value") %>% 
  mutate(cngr = factor(cngr, levels = c("Z_ms", "Z_m6s", "Z_m6m", "Z_m6ms"))) %>% 
  mutate(cngr = recode_factor(cngr, "Z_ms" = "1", "Z_m6s" = "2", "Z_m6m" = "3", "Z_m6ms" = "4")) %>% # trick: numeric to plot geom_line
  mutate(cngr = as.numeric(cngr)) %>% 
  ggplot(aes(x = cngr, y = value, color = ID)) + 
  geom_line(alpha = 0.2) + geom_point(alpha = 0.6) + 
  scale_color_manual(values = rep(pal[2], 212)) + 
  scale_x_continuous(labels = c("Mod-Sct", "Subg-Sct", "Subg-Mod", "Subg-\nMod-Sct")) + 
  xlab("") + ylab("") + facet_wrap(~Code, nrow = 4, scales = "free_y") + theme_cowplot() + 
  theme(legend.position = "none", plot.margin = unit(c(0,0.5,0,0), "cm")); M

miss %>% filter(Code %in% c("Pol-P", "SD-P", "F-Ane", "Ant-P")) %>% 
  filter(zscore %in% c("Z_ms", "Z_m6s", "Z_m6m", "Z_m6ms")) %>% 
  filter(value == T) %>% ungroup() %>% select(n) %>% sum()

# ----- NORMALISED ZSCORE ----- #

nA <- zdata %>% filter(Code %in% c("Pa-Ho", "Ph-B", "Her-P", "Pr-Pr")) %>%
  select(ID, Code, contains("Norm")) %>%
  pivot_longer(3:6, names_to = "cngr", values_to = "value") %>% 
  mutate(cngr = factor(cngr, levels = c("NormZ_ms", "NormZ_m6s", "NormZ_m6m", "NormZ_m6ms"))) %>% 
  mutate(cngr = recode_factor(cngr, "NormZ_ms" = "1", "NormZ_m6s" = "2", "NormZ_m6m" = "3", "NormZ_m6ms" = "4")) %>% # trick: numeric to plot geom_line
  mutate(cngr = as.numeric(cngr)) %>% 
  ggplot(aes(x = cngr, y = value, color = ID)) + 
  geom_line(alpha = 0.2) + geom_point(alpha = 0.6) + 
  scale_color_manual(values = rep(pal[1], 154)) + 
  scale_x_continuous(labels = c("Mod-Sct", "Subg-Sct", "Subg-Mod", "Subg-\nMod-Sct")) + 
  xlab("") + ylab("Congruence coefficients (normalised z-scores)") +
  facet_wrap(~Code, nrow = 4, scales = "free_y") + theme_cowplot() + 
  theme(legend.position = "none", plot.margin = unit(c(0,0.5,0,0), "cm")); nA

nM <- zdata %>% filter(Code %in% c("Pol-P", "SD-P", "F-Ane", "Ant-P")) %>%
  select(ID, Code, contains("Norm")) %>%
  pivot_longer(3:6, names_to = "cngr", values_to = "value") %>% 
  mutate(cngr = factor(cngr, levels = c("NormZ_ms", "NormZ_m6s", "NormZ_m6m", "NormZ_m6ms"))) %>% 
  mutate(cngr = recode_factor(cngr, "NormZ_ms" = "1", "NormZ_m6s" = "2", "NormZ_m6m" = "3", "NormZ_m6ms" = "4")) %>% # trick: numeric to plot geom_line
  mutate(cngr = as.numeric(cngr)) %>% 
  ggplot(aes(x = cngr, y = value, color = ID)) + 
  geom_line(alpha = 0.2) + geom_point(alpha = 0.6) + 
  scale_color_manual(values = rep(pal[2], 212)) + 
  scale_x_continuous(labels = c("Mod-Sct", "Subg-Sct", "Subg-Mod", "Subg-\nMod-Sct")) + 
  xlab("") + ylab("") + facet_wrap(~Code, nrow = 4, scales = "free_y") + theme_cowplot() + 
  theme(legend.position = "none", plot.margin = unit(c(0,0.5,0,0), "cm")); nM

csp <- plot_grid(A, M, nA, nM, nrow = 2)
title <- ggdraw() + draw_label("Congruence Significance Profile", fontface = 'bold')
plot_grid(title, csp, ncol = 1, rel_heights = c(0.05, 1, 1))
