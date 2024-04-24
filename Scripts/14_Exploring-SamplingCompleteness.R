# -------------- Manuscript figures: coevolutionary effects ---------------

# Author: Kate P Maia
# Checked: 04/2024

# Explores main structural and dynamics results on cropped dataset: only including networks which attempted to quantify sampling completeness. 

# --------------------- Loading library, code and data --------------------

setwd("C:/Users/Kate Maia/Documents/Maia-Guimaraes_HierarchyCoevoUnits")

library(tidyverse)
library(RColorBrewer)
library(cowplot)
#source("./Scripts/functions/zscore.R")

sampling <- read.table("../../Dropbox/Kate_Manuscripts/Hierarchical Structure/EL_Submission2/TableS1_Sampling.txt", header = TRUE)

net_struct <- read.table("./Outputs/05_Net-Level_Struct.txt", header = T, sep = "\t")

part_SM <- read.table("./Outputs/09_Partition-Effects-SctMd-Dataframe.txt", header = T, sep = "\t")
part_MM <- read.table("./Outputs/10_Partition-Effects-MdMt-Dataframe.txt", header = T, sep = "\t")

# Combining effects dataframes (partSM and partMM) - CHECKED
partition_df <- part_SM %>% filter(effmeasure != "teff") %>% 
  left_join(part_MM, by = c("ID", "IntType", "effmeasure", "m"))
head(partition_df); nrow(part_SM %>% filter(effmeasure != "teff")) == nrow(partition_df)

# Adding interaction type and sign information to all datasets
intlevels <- c("Host-Parasite", "Phage-Bacteria", "Plant-Herbivore", "Food Webs", "Pollination", "Seed Dispersal", "Anemone-Fish", "Plant-Ant")

net_struct$IntSign <- net_struct$Code <- net_struct$IntType <- factor(net_struct$IntType, levels = intlevels)
levels(net_struct$IntType) <- c("Parasite-Host", "Phage-Bacteria", "Herbivore-Plant", "Predator-Prey", "Pollinator-Plant", "Seed Disperser-Plant", "Fish-Anemone", "Ant-Plant")
levels(net_struct$Code) <- c("Pa-Ho", "Ph-B", "Her-P", "Pr-Pr", "Pol-P", "SD-P", "F-Ane", "Ant-P")
levels(net_struct$IntSign) <- rep(c("A", "M"), each = 4)

partition_df$Code <- partition_df$IntType <- factor(partition_df$IntType, levels = intlevels)
levels(partition_df$IntType) <- levels(net_struct$IntType)
levels(partition_df$Code) <- levels(net_struct$Code)
  
# multLC & starLC (to be removed plot by plot)
multLC <- c("Herb_CuevasReyes", "IZ_PA_DIM100", "IZ_PA_PAcamp", "M_AF_002_02", "M_AF_002_11")
starLC <- c("IZ_PA_DIM10", "IZ_PA_PA10", "M_AF_002_07", "M_AF_002_16", "M_PL_061_33")

strpal <- c("#D85C00", "#31489F")
effpal <- rev(c(brewer.pal(8, "Set2")[c(3:5)], "grey90"))

# -------------------------------------------------------------------------
# -------------------------- Summary of sampling --------------------------

sampling %>% filter(Sampling == "Yes") %>% group_by(IntType) %>% count()
sampled <- sampling %>% filter(Sampling == "Yes") %>% select(ID); sampled <- sampled$ID

net_struct %>% filter(ID %in% sampled) %>% # counting number of components
  group_by(IntType, NComp) %>% count() %>% arrange(NComp)

net_struct %>% filter(ID %in% sampled) %>% # size of largest components
  mutate(LargeCompSize = LargeCompSize/(NRow + NCol)) %>% 
  group_by(IntType) %>% summarise(mean = mean(LargeCompSize), sd = sd(LargeCompSize))

sampled[sampled %in% c(multLC)] # 2: 2 Ant-P
sampled[sampled %in% c(starLC)] # 3: 1 Pol-P, 2 Ant-P

# -------------------------------------------------------------------------
# -------------------------- Exploring STRUCTURE ---------------------------

data <- net_struct %>% filter(!ID %in% c(multLC, starLC)) # 366 ID
data <- data %>% select(1,35,36,15:31); head(data); dim(data) # facilitate visualization

# -------------------------------------------------------------------------
# --------------- Congruence: zscore and normalised zscore ----------------

# module-sector congruence
msc <- data %>% select(ID, Code, MS_congr, MS_mean, MS_sd); head(msc)
msc <- zscore(msc); msc %>% filter(is.na(Z) | is.na(NormZ))

msdata <- data %>% select(ID, IntSign, Code, contains("MS_")) # saves full mod-sect data (N = 366) before cropping 
data <- data %>% filter(n_m6 != 0); dim(data) # if no subgraphs, no point in testing their hierarchy 

# subgraph-sector congruence
m6sc <- data %>% select(ID, Code, n_m6_sect, mean_m6_sect, sd_m6_sect); head(m6sc)
m6sc <- zscore(m6sc); m6sc %>% filter(is.na(Z)) %>% count() # 6 NaN (0/0)
dim(m6sc) # 0 removed due to Inf: Sd 0 become Inf: 0/0 is Nan; 1/0 is Inf

# subgraph-module congruence
m6mc <- data %>% select(ID, Code, n_m6_mod, mean_m6_mod, sd_m6_mod); head(m6mc)
m6mc %>% filter(is.na(n_m6_mod) | is.na(mean_m6_mod) | is.na(sd_m6_mod))
m6mc <- zscore(m6mc); m6mc %>% filter(is.na(Z)) # 24 Sd = 0
dim(m6mc) # 0 removed due to Inf

# subgraph-module-sector congruence ("double hierarchy")
m6msc <- data %>% filter(ID != "M_PL_062") %>% # additionally removes M_PL_062 - too large
  select(ID, Code, n_m6_ms, mean_m6_ms, sd_m6_ms); head(m6msc); dim(m6msc) # 344
m6msc %>% filter(is.na(n_m6_ms) | is.na(mean_m6_ms) | is.na(sd_m6_ms))
m6msc <- zscore(m6msc); m6msc %>% filter(is.na(Z)) # 36 Sd = 0
dim(m6msc) # 2 removed due to Inf

zdata <- left_join(msc, m6sc, by = c("ID", "Code"), suffix = c("_ms", "_m6s")) %>% 
  left_join(m6mc, by = c("ID", "Code")) %>% 
  left_join(m6msc, by = c("ID", "Code"), suffix = c("_m6m", "_m6ms")); dim(zdata)

# -------------------------------------------------------------------------
# ---------- Fig SX: Plotting module-sector congruence results ------------

msdata <- msdata %>% filter(ID %in% sampled); dim(msdata) # 79 (84 sampled - 5 starLD/multLC)

msdata %>% summarise(mean = mean(MS_congr), sd = sd(MS_congr), min = min(MS_congr))

pa <- msdata %>% select(ID, IntSign, Code, contains("MS_")) %>% # N = 79
  ggplot(aes(x = Code, y = MS_congr, fill = IntSign)) +
  geom_boxplot(color = c(rep(strpal[1], 1), rep(strpal[2], 3)), alpha = 0.5) +
  xlab("") + ylab("Module-sector congruence") +
  coord_cartesian(ylim = c(0.50, 1.00)) + scale_fill_manual(values = strpal) + 
  theme_classic() + theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 14), legend.position = "none"); pa

pb <- msdata %>% select(ID, IntSign, Code, MS_pval) %>% # N = 79 
  mutate(plogic = ifelse(MS_pval <= 0.05, TRUE, FALSE)) %>% 
  mutate(plogic = factor(plogic, levels = c(FALSE, TRUE))) %>% 
  group_by(Code, plogic) %>% summarise(count = n()) %>% 
  ggplot(aes(x = Code, y = count, fill = Code, alpha = plogic)) + 
  geom_bar(position = "fill", stat = "identity") + # position = "stack"
  xlab("") + ylab("Proportion of networks") + scale_alpha_manual(values = c(0.2, 0.9)) + 
  scale_fill_manual(values = c(rep(strpal[1], 1), rep(strpal[2], 3))) +
  ggtitle("Module-sector congruence: null model") +
  theme_classic() + theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 14), plot.title = element_text(size = 14), legend.position = "none")
pb <- pb + annotate("text", x = 1:4, y = 1.02, label = as.numeric(table(msdata$Code)[table(msdata$Code) != 0])); pb

pc <- zdata %>% filter(ID %in% sampled) %>% # 79
  ggplot() + geom_histogram(aes(x = Z_ms), alpha = 0.6) + 
  xlab("Module-sector congruence (z-score)") + ylab("Frequecy") + theme_cowplot(); pc

cowplot::plot_grid(pa, pb, pc, labels = c("a)", "b)", "c)"), nrow = 1, label_size = 12, label_fontfamily = "", label_fontface = "plain")
#ggsave("../../Dropbox/Kate_Manuscripts/Hierarchical Structure/EL_Submission2/FigSI_Sampling_MSCongr.svg", width = 33, height = 11, units = "cm", bg = "white")

# -------------------------------------------------------------------------
# ------------ Fig SX: Plotting subgraphs in modules/sectors -------------- 

net_struct %>% filter(!ID %in% c(starLC, multLC), ID %in% sampled) %>% # N = 79
  mutate(n_m6 = n_m6 == 0) %>% group_by(IntType, n_m6) %>% count() %>% arrange(n_m6)

data %>% filter(ID %in% sampled) %>% # N = 74 with n_m6
  select(ID, IntSign, Code, n_m6, n_m6_sect) %>% mutate(In = n_m6_sect/n_m6) %>% 
  summarise(mean = mean(In), sd = sd(In), min = min(In), max = max(In))
data %>% filter(ID %in% sampled) %>% # N = 74 with n_m6
  select(ID, IntSign, Code, p_m6_sect) %>% mutate(p = p_m6_sect <= 0.05) %>% count(p)

pa <- data %>% filter(ID %in% sampled) %>% # N = 74 with n_m6
  select(ID, IntSign, Code, n_m6, n_m6_sect) %>% 
  mutate(In = n_m6_sect/n_m6) %>%
  ggplot(aes(x = Code, y = In, col = IntSign, fill = IntSign)) + geom_boxplot(alpha = 0.5) +
  scale_color_manual(values = strpal) + scale_fill_manual(values = strpal) + 
  xlab("") + ylab("Subgraph-sector congruence") +
  theme_classic() + theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 14), legend.position = "none"); pa

label <- data %>% filter(ID %in% sampled) %>% group_by(Code) %>% count()
pb <- data %>% filter(ID %in% sampled) %>% # N = 74 with n_m6
  select(ID, IntSign, Code, p_m6_sect) %>%
  mutate(p_m6_sect = ifelse(p_m6_sect <= 0.05, TRUE, FALSE)) %>% 
  mutate(p_m6_sect = factor(p_m6_sect, levels = c(FALSE, TRUE))) %>% 
  group_by(Code, p_m6_sect) %>% summarise(count = n()) %>%
  ggplot(aes(x = Code, y = count, fill = Code, alpha = p_m6_sect)) + 
  geom_bar(position = "fill", stat = "identity") + # position = "stack"
  xlab("") + ylab("Proportion of networks") + ggtitle("Subgraph-sector congruence: null model") +
  scale_alpha_manual(values = c(0.2, 0.9)) + scale_fill_manual(values = c(strpal[1], rep(strpal[2], 3))) +
  theme_classic() + theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 14), plot.title = element_text(size = 14), legend.position = "none")
pb <- pb + annotate("text", x = 1:4, y = 1.02, label = label$n); pb

pc <- zdata %>% filter(ID %in% sampled) %>% # N = 79 
  filter(!is.na(Z_m6s)) %>% # N = 74 (see above) 
  ggplot() + geom_histogram(aes(x = Z_m6s), alpha = 0.6) +
  xlab("Subgraph-sector congruence (z-score)") + ylab("Frequecy") + 
  theme_cowplot() + theme(axis.title.x = element_text(size = 12)); pc

data %>% filter(ID %in% sampled) %>% # N = 74 with n_m6
  select(ID, IntSign, Code, n_m6, n_m6_mod) %>% mutate(In = n_m6_mod/n_m6) %>% 
  summarise(mean = mean(In), sd = sd(In), min = min(In), max = max(In))
data %>% filter(ID %in% sampled) %>% # N = 74 with n_m6
  select(ID, IntSign, Code, p_m6_mod) %>% mutate(p = p_m6_mod <= 0.05) %>% count(p)

pd <- data %>% filter(ID %in% sampled) %>% # N = 74 with n_m6
  select(ID, IntSign, Code, n_m6, n_m6_mod) %>%
  mutate(In = n_m6_mod/n_m6) %>%
  ggplot(aes(x = Code, y = In, col = IntSign, fill = IntSign)) + geom_boxplot(alpha = 0.5) +
  scale_color_manual(values = strpal) + scale_fill_manual(values = strpal) + 
  xlab("") + ylab("Subgraph-module congruence") + theme_classic() + 
  theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 14), legend.position = "none"); pd

pe <- data %>% filter(ID %in% sampled) %>%  # N = 74 with n_m6
  select(ID, IntSign, Code, p_m6_mod) %>%
  mutate(p_m6_mod = ifelse(p_m6_mod <= 0.05, TRUE, FALSE)) %>% 
  mutate(p_m6_mod = factor(p_m6_mod, levels = c(FALSE, TRUE))) %>% 
  group_by(Code, p_m6_mod) %>% summarise(count = n()) %>% 
  ggplot(aes(x = Code, y = count, fill = Code, alpha = p_m6_mod)) + 
  geom_bar(position = "fill", stat = "identity") + # position = "stack"
  xlab("") + ylab("Proportion of networks") + ggtitle("Subgraph-module congruence: null model") +
  scale_alpha_manual(values = c(0.2, 0.9)) + scale_fill_manual(values = c(strpal[1], rep(strpal[2], 3))) +
  theme_classic() + theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 14), plot.title = element_text(size = 14), legend.position = "none"); pe

pf <- zdata %>% filter(ID %in% sampled) %>% # N = 79
  filter(!is.na(Z_m6m)) %>% # N = 69 (see above) 
  ggplot() + geom_histogram(aes(x = Z_m6m), alpha = 0.6) + 
  xlab("Subgraph-module congruence (z-score)") + ylab("Frequecy") + 
  theme_cowplot() + theme(axis.title.x = element_text(size = 12)); pf

cowplot::plot_grid(pa, pb, pc, pd, pe, pf, labels = c("a)", "b)", "c)", "d)", "e)", "f)"), nrow = 2, label_size = 14, label_fontfamily = "", label_fontface = "plain")
#ggsave("../../Dropbox/Kate_Manuscripts/Hierarchical Structure/EL_Submission2/FigSI_Sampling_SubgrCongr.svg", width = 33, height = 22, units = "cm", bg = "white")

# -------------------------------------------------------------------------
# --------------- Fig SX: Subgraph-Module-Sector Congruence ---------------

data %>% filter(ID %in% sampled) %>% # N = 74 with n_m6
  select(ID, IntSign, Code, n_m6, n_m6_ms) %>% mutate(In = n_m6_ms/n_m6) %>% 
  summarise(mean = mean(In), sd = sd(In), min = min(In), max = max(In))
data %>% filter(ID %in% sampled) %>% # N = 74 with n_m6
  select(ID, IntSign, Code, p_m6_ms) %>% mutate(p = p_m6_ms <= 0.05) %>% count(p)

pa <- data %>% filter(ID %in% sampled) %>% # N = 74 with n_m6
  select(Code, IntSign, ID, n_m6, n_m6_ms) %>% # all network types
  mutate(In = n_m6_ms/n_m6) %>%
  ggplot(aes(x = Code, y = In, col = IntSign, fill = IntSign)) + geom_boxplot(alpha = 0.5) +
  scale_color_manual(values = strpal) + scale_fill_manual(values = strpal) +
  xlab("") + ylab("Subgraphs-module-sector congruence") + theme_classic() +
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 12), legend.position = "none"); pa

label <- data %>% filter(ID %in% sampled) %>% group_by(Code) %>% count()
pb <- data %>% filter(ID %in% sampled) %>% # N = 74 with n_m6
  mutate(p_m6_ms = ifelse(p_m6_ms <= 0.05, TRUE, FALSE)) %>% 
  mutate(p_m6_ms = factor(p_m6_ms, levels = c(FALSE, TRUE))) %>% 
  group_by(Code, p_m6_ms) %>% summarise(count = n()) %>% 
  ggplot(aes(x = Code, y = count, fill = Code, alpha = p_m6_ms)) + 
  geom_bar(position = "fill", stat = "identity") + # position = "stack"
  xlab("") + ylab("Proportion of networks") +
  scale_alpha_manual(values = c(0.2, 0.9)) + scale_fill_manual(values = c(strpal[1], rep(strpal[2], 3))) +
  ggtitle("Subgraphs-module-sector congruence: null model") + theme_classic() +
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 12), plot.title = element_text(size = 14), legend.position = "none")
pb <- pb + annotate("text", x = 1:4, y = 1.02, label = label$n); pb

cowplot::plot_grid(pa, pb, labels = c("a)", "b)"), ncol = 2, label_size = 14, label_fontfamily = "", label_fontface = "plain")
#ggsave("../../Dropbox/Kate_Manuscripts/Hierarchical Structure/EL_Submission2/FigSI_Sampling_DoubleHier.svg.svg", width = 26, height = 13, units = "cm", bg = "white")

# -------------------------------------------------------------------------
# -------------------------- Exploring DYNAMICS ---------------------------

partition_df <- partition_df %>% filter(ID %in% sampled) # 79 (84 sampled - 5 starLD/multLC)

# -------------------------------------------------------------------------
# -------- Fig SX: Plotting module-subgraph partition of effects ----------

partition_df %>% select(ID, Code, effmeasure, m, inMdinMt, inMdbetMt, betMdinMt, betMdbetMt) %>%
  mutate(tot = inMdinMt + inMdbetMt + betMdinMt + betMdbetMt) %>% 
  filter(!is.na(tot)) %>% # removes 2 large ID (Herb_Pearse&Alterm, M_PL_062)
  mutate(inMdinMt = inMdinMt/tot, inMdbetMt = inMdbetMt/tot, betMdinMt = betMdinMt/tot, betMdbetMt = betMdbetMt/tot) %>% 
  group_by(Code, effmeasure, m) %>% # 8Code*2effmeasure*3m = 48groups
  summarise(inMdinMt = mean(inMdinMt), inMdbetMt = mean(inMdbetMt), betMdinMt = mean(betMdinMt), betMdbetMt = mean(betMdbetMt)) %>% 
  mutate(effmeasure = recode(effmeasure, deff = "Direct", ieff = "Indirect")) %>% 
  pivot_longer(cols = 4:7, names_to = "Partition", values_to = "value") %>% 
  mutate(Partition = factor(Partition, levels = c("betMdbetMt", "betMdinMt", "inMdbetMt", "inMdinMt"))) %>%
  ggplot(aes(x = Code, y = value, fill = Partition)) + geom_bar(stat = "identity") +
  labs(x = "", y = "Proportion of effects") + theme_classic() +
  theme(legend.position = "none") + scale_fill_manual(values = effpal) +
  facet_grid(rows = vars(m), cols = vars(effmeasure))
#ggsave("../../Dropbox/Kate_Manuscripts/Hierarchical Structure/EL_Submission2/FigSI_Sampling_PartitionModSubg.svg", width = 6.2, height = 6.2)

# -------------------------------------------------------------------------
# --------- Fig SX: Plotting sector-module partition of effects -----------

partition_df %>% select(ID, Code, effmeasure, m, inSinM, inSbetM, betSinM, betSbetM) %>%
  mutate(tot = inSinM + inSbetM + betSinM + betSbetM) %>%
  mutate(inSinM = inSinM/tot, inSbetM = inSbetM/tot, betSinM = betSinM/tot, betSbetM = betSbetM/tot) %>%
  group_by(Code, effmeasure, m) %>% # 8Code*2effmeasure*3m = 48groups
  summarise(inSinM = mean(inSinM), inSbetM = mean(inSbetM), betSinM = mean(betSinM), betSbetM = mean(betSbetM)) %>% 
  mutate(effmeasure = recode(effmeasure, deff = "Direct", ieff = "Indirect")) %>% 
  pivot_longer(cols = 4:7, names_to = "Partition", values_to = "value") %>% 
  mutate(Partition = factor(Partition, levels = c("betSbetM", "betSinM", "inSbetM", "inSinM"))) %>%
  ggplot(aes(x = Code, y = value, fill = Partition)) + geom_bar(stat = "identity") +
  labs(x = "", y = "Proportion of effects") + theme_classic() +
  theme(legend.position = "none") + scale_fill_manual(values = effpal) +
  facet_grid(rows = vars(m), cols = vars(effmeasure))
#ggsave("../../Dropbox/Kate_Manuscripts/Hierarchical Structure/EL_Submission2/FigSI_Sampling_PartitionSectMod.svg", width = 6.2, height = 6.2)

# -------------------------------------------------------------------------
# --------- Results on direct effects within groups across scales ---------

deffdf <- partition_df %>% filter(effmeasure == "deff", m == 0.1) %>% # deff: indif to m when prop
  mutate(totMM = inMdinMt + inMdbetMt + betMdinMt + betMdbetMt,
         totSM = inSinM + inSbetM + betSinM + betSbetM) %>%
  mutate(inMdinMt = inMdinMt/totMM, inMdbetMt = inMdbetMt/totMM, betMdinMt = betMdinMt/totMM, betMdbetMt = betMdbetMt/totMM,
         inSinM = inSinM/totSM, inSbetM = inSbetM/totSM, betSinM = betSinM/totSM, betSbetM = betSbetM/totSM) %>%
  select(-totMM, -totSM) %>% 
  mutate(inMt = inMdinMt + betMdinMt, inMd = inMdinMt + inMdbetMt, inM = inSinM + betSinM, inS = inSinM + inSbetM) 

deffdf %>% select(inMt, inMd, inM, inS) %>%
  summarise(meinMt = mean(inMt, na.rm = TRUE)*100, sdinMt = sd(inMt, na.rm = TRUE)*100, 
            meinMd = mean(inMd, na.rm = TRUE)*100, sdinMd = sd(inMd, na.rm = TRUE)*100,
            meinM = mean(inM)*100, sdinM = sd(inM)*100, 
            meinS = mean(inS)*100, sdinS = sd(inS)*100)

# -------------------------------------------------------------------------
# -- Results on (in)direct effects and hierarchical structure (Subg/Mod) --

partition_df %>% # MEAN
  mutate(totMM = inMdinMt + inMdbetMt + betMdinMt + betMdbetMt) %>%
  mutate(inMdinMt = inMdinMt/totMM, inMdbetMt = inMdbetMt/totMM, betMdinMt = betMdinMt/totMM, betMdbetMt = betMdbetMt/totMM) %>%
  select(-totMM, -Code) %>% group_by(effmeasure, m) %>%
  summarise(inMdinMt = mean(inMdinMt, na.rm = TRUE), inMdbetMt = mean(inMdbetMt, na.rm = TRUE), 
            betMdinMt = mean(betMdinMt, na.rm = TRUE), betMdbetMt = mean(betMdbetMt, na.rm = TRUE)) %>% 
  pivot_wider(names_from = "effmeasure", values_from = 3:6) %>% 
  mutate(inMdinMt_D = inMdinMt_deff*100, inMdbetMt_D = inMdbetMt_deff*100, 
         betMdinMt_D = betMdinMt_deff*100, betMdbetMt_D = betMdbetMt_deff*100, 
         inMdinMt_I = inMdinMt_ieff*100, inMdbetMt_I = inMdbetMt_ieff*100,  
         betMdinMt_I = betMdinMt_ieff*100, betMdbetMt_I = betMdbetMt_ieff*100) %>% 
  select(-contains("deff"), -contains("ieff"))

partition_df %>% # SD
  mutate(totMM = inMdinMt + inMdbetMt + betMdinMt + betMdbetMt) %>%
  mutate(inMdinMt = inMdinMt/totMM, inMdbetMt = inMdbetMt/totMM, betMdinMt = betMdinMt/totMM, betMdbetMt = betMdbetMt/totMM) %>%
  select(-totMM, -Code) %>% group_by(effmeasure, m) %>%
  summarise(inMdinMt = sd(inMdinMt, na.rm = TRUE), inMdbetMt = sd(inMdbetMt, na.rm = TRUE), 
            betMdinMt = sd(betMdinMt, na.rm = TRUE), betMdbetMt = sd(betMdbetMt, na.rm = TRUE)) %>% 
  pivot_wider(names_from = "effmeasure", values_from = 3:6) %>% 
  mutate(inMdinMt_D = inMdinMt_deff*100, inMdbetMt_D = inMdbetMt_deff*100, 
         betMdinMt_D = betMdinMt_deff*100, betMdbetMt_D = betMdbetMt_deff*100, 
         inMdinMt_I = inMdinMt_ieff*100, inMdbetMt_I = inMdbetMt_ieff*100,  
         betMdinMt_I = betMdinMt_ieff*100, betMdbetMt_I = betMdbetMt_ieff*100) %>% 
  select(-contains("deff"), -contains("ieff"))

# -------------------------------------------------------------------------
# -- Results on (in)direct effects and hierarchical structure (Mod/Sect) --

partition_df %>% # MEAN
  mutate(totSM = inSinM + inSbetM + betSinM + betSbetM) %>%
  mutate(inSinM = inSinM/totSM, inSbetM = inSbetM/totSM, betSinM = betSinM/totSM, betSbetM = betSbetM/totSM) %>%
  select(-totSM, -Code) %>% group_by(effmeasure, m) %>%
  summarise(inSinM = mean(inSinM), inSbetM = mean(inSbetM), betSinM = mean(betSinM), betSbetM = mean(betSbetM)) %>% 
  pivot_wider(names_from = "effmeasure", values_from = 3:6) %>% 
  mutate(inSinM_D = inSinM_deff*100, inSbetM_D = inSbetM_deff*100, 
         betSinM_D = betSinM_deff*100, betSbetM_D = betSbetM_deff*100, 
         inSinM_I = inSinM_ieff*100, inSbetM_I = inSbetM_ieff*100,  
         betSinM_I = betSinM_ieff*100, betSbetM_I = betSbetM_ieff*100) %>% 
  select(-contains("deff"), -contains("ieff"))

partition_df %>% # SD
  mutate(totSM = inSinM + inSbetM + betSinM + betSbetM) %>%
  mutate(inSinM = inSinM/totSM, inSbetM = inSbetM/totSM, betSinM = betSinM/totSM, betSbetM = betSbetM/totSM) %>%
  select(-totSM, -Code) %>% group_by(effmeasure, m) %>%
  summarise(inSinM = sd(inSinM), inSbetM = sd(inSbetM), betSinM = sd(betSinM), betSbetM = sd(betSbetM)) %>% 
  pivot_wider(names_from = "effmeasure", values_from = 3:6) %>% 
  mutate(inSinM_D = inSinM_deff*100, inSbetM_D = inSbetM_deff*100, 
         betSinM_D = betSinM_deff*100, betSbetM_D = betSbetM_deff*100, 
         inSinM_I = inSinM_ieff*100, inSbetM_I = inSbetM_ieff*100,  
         betSinM_I = betSinM_ieff*100, betSbetM_I = betSbetM_ieff*100) %>% 
  select(-contains("deff"), -contains("ieff"))
