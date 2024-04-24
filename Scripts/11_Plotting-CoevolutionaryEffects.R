# -------------- Manuscript figures: coevolutionary effects ---------------

# Author: Kate P Maia
# Checked: 12/2023

# Reads data generated in previous scripts to plot manuscript figures on coevolutionary effects across groups and the hierarchical structure of ecological networks. 

# --------------------- Loading library, code and data --------------------

setwd("C:/Users/Kate Maia/Documents/Maia-Guimaraes_HierarchyCoevoUnits")

library(tidyverse)
library(reshape2)
library(cowplot)
library(RColorBrewer)
library(GGally)
library(grid)
library(gridExtra)
source("./Scripts/Functions/lines_from_points_grp.R")

net_struct <- read.table("./Outputs/05_Net-Level_Struct.txt", header = T, sep = "\t")

part_SM <- read.table("./Outputs/09_Partition-Effects-SctMd-Dataframe.txt", header = T, sep = "\t")
part_MM <- read.table("./Outputs/10_Partition-Effects-MdMt-Dataframe.txt", header = T, sep = "\t")

# Combining effects dataframes (partSM and partMM) - CHECKED
partition_df <- part_SM %>% filter(effmeasure != "teff") %>% 
  left_join(part_MM, by = c("ID", "IntType", "effmeasure", "m"))
head(partition_df); nrow(part_SM %>% filter(effmeasure != "teff")) == nrow(partition_df)

intlevels <- c("Host-Parasite", "Phage-Bacteria", "Plant-Herbivore", "Food Webs", "Pollination", "Seed Dispersal", "Anemone-Fish", "Plant-Ant")

partition_df$IntType <- factor(partition_df$IntType, levels = intlevels)
net_struct$IntType <- factor(net_struct$IntType, levels = intlevels)
levels(net_struct$IntType) <- levels(partition_df$IntType) <- c("Parasite-Host", "Phage-Bacteria", "Herbivore-Plant", "Predator-Prey", "Pollinator-Plant", "Seed Disperser-Plant", "Fish-Anemone", "Ant-Plant")

partition_df$Code <- partition_df$IntType; net_struct$Code <- net_struct$IntType
levels(net_struct$Code) <- levels(partition_df$Code) <- c("Pa-Ho", "Ph-B", "Her-P", "Pr-Pr", "Pol-P", "SD-P", "F-Ane", "Ant-P")

table(partition_df$IntType, partition_df$m, partition_df$effmeasure)

pal <- rev(c(brewer.pal(8, "Set2")[c(3:5)], "grey90"))

# -------------------------------------------------------------------------
# --------- Fig 4: Plotting module-subgraph partition of effects ----------

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
  theme(legend.position = "none") + scale_fill_manual(values = pal) +
  facet_grid(rows = vars(m), cols = vars(effmeasure))

# -------------------------------------------------------------------------
# ---------- Fig 5: Plotting sector-module partition of effects -----------

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
  theme(legend.position = "none") + scale_fill_manual(values = pal) +
  facet_grid(rows = vars(m), cols = vars(effmeasure))

# -------------------------------------------------------------------------
# --- Fig X: Module-subgraph partition effect as significance profiles ----

spdata <- partition_df %>% select(ID, Code, effmeasure, m, contains("Mt")) %>%
  mutate(tot = inMdinMt + inMdbetMt + betMdinMt + betMdbetMt) %>% 
  filter(!is.na(tot)) %>% # removes 2 large ID (Herb_Pearse&Alterm, M_PL_062)
  mutate(inMdinMt = inMdinMt/tot, inMdbetMt = inMdbetMt/tot, betMdinMt = betMdinMt/tot, betMdbetMt = betMdbetMt/tot) %>% select(-tot) %>% 
  mutate(effmeasure = recode(effmeasure, deff = "Direct", ieff = "Indirect")) %>%
  pivot_longer(cols = 5:8, names_to = "Partition", values_to = "value") %>% 
  mutate(Partition = factor(Partition, levels = c("inMdinMt", "inMdbetMt", "betMdinMt", "betMdbetMt"))) %>%
  mutate(Partition = recode_factor(Partition, "inMdinMt" = "1", "inMdbetMt" = "2", "betMdinMt" = "3", "betMdbetMt" = "4")) %>% # trick: numeric to plot geom_line
  mutate(Partition = as.numeric(Partition))

lab <- c( "inSg\ninM", "betSg\ninM", "inSg\nbetM", "betSg\nbetM")

paho <- spdata %>% filter(Code == "Pa-Ho") %>% # 1944 = 4 * 6 * 81ID
  ggplot(aes(x = Partition, y = value, color = ID)) + 
  geom_line(alpha = 0.2) + geom_point(alpha = 0.2) + xlab("") + ylab("") + 
  scale_color_manual(values = rep("#d94801", 81)) + 
  scale_x_continuous(labels = lab) + scale_y_continuous(limits = c(0, 1)) + 
  theme_classic() + theme(legend.position = "none") + ggtitle("Parasite-Host") + 
  facet_grid(rows = vars(m), cols = vars(effmeasure)); paho
  
phb <- spdata %>% filter(Code == "Ph-B") %>% # 912 = 4 * 6 * 38ID
  ggplot(aes(x = Partition, y = value, color = ID)) + 
  geom_line(alpha = 0.2) + geom_point(alpha = 0.2) + xlab("") + ylab("") + 
  scale_color_manual(values = rep("#d94801", 38)) + 
  scale_x_continuous(labels = lab) + scale_y_continuous(limits = c(0, 1)) + 
  theme_classic() + theme(legend.position = "none") + ggtitle("Phage-Bacteria") +
  facet_grid(rows = vars(m), cols = vars(effmeasure)); phb
  
herp <- spdata %>% filter(Code == "Her-P") %>% # 744 = 4 * 6 * 31ID
  ggplot(aes(x = Partition, y = value, color = ID)) + 
  geom_line(alpha = 0.2) + geom_point(alpha = 0.2) + xlab("") + ylab("") + 
  scale_color_manual(values = rep("#d94801", 32)) + 
  scale_x_continuous(labels = lab) + scale_y_continuous(limits = c(0, 1)) + 
  theme_classic() + theme(legend.position = "none") + ggtitle("Herbivore-Plant") +
  facet_grid(rows = vars(m), cols = vars(effmeasure)); herp

prpr <- spdata %>% filter(Code == "Pr-Pr") %>% # 72 = 4 * 6 * 3ID
  ggplot(aes(x = Partition, y = value, color = ID)) + 
  geom_line(alpha = 0.2) + geom_point(alpha = 0.2) + xlab("") + ylab("") + 
  scale_color_manual(values = rep("#d94801", 3)) + 
  scale_x_continuous(labels = lab) + scale_y_continuous(limits = c(0, 1)) + 
  theme_classic() + theme(legend.position = "none") + ggtitle("Predator-Prey") +
  facet_grid(rows = vars(m), cols = vars(effmeasure)); prpr

polp <- spdata %>% filter(Code == "Pol-P") %>% # 3648 = 4 * 6 * 152ID
  ggplot(aes(x = Partition, y = value, color = ID)) + 
  geom_line(alpha = 0.2) + geom_point(alpha = 0.2) + xlab("") + ylab("") + 
  scale_color_manual(values = rep("#2171b5", 152)) + 
  scale_x_continuous(labels = lab) + scale_y_continuous(limits = c(0, 1)) + 
  theme_classic() + theme(legend.position = "none") + ggtitle("Pollinator-Plant") +
  facet_grid(rows = vars(m), cols = vars(effmeasure)); polp

sdp <- spdata %>% filter(Code == "SD-P") %>% # 936 = 4 * 6 * 39ID
  ggplot(aes(x = Partition, y = value, color = ID)) + 
  geom_line(alpha = 0.2) + geom_point(alpha = 0.2) + xlab("") + ylab("") + 
  scale_color_manual(values = rep("#2171b5", 39)) + 
  scale_x_continuous(labels = lab) + scale_y_continuous(limits = c(0, 1)) + 
  theme_classic() + theme(legend.position = "none") + ggtitle("Seed Disperser-Plant") +
  facet_grid(rows = vars(m), cols = vars(effmeasure)); sdp

fane <- spdata %>% filter(Code == "F-Ane") %>% # 312 = 4 * 6 * 13ID
  ggplot(aes(x = Partition, y = value, color = ID)) + 
  geom_line(alpha = 0.2) + geom_point(alpha = 0.2) + xlab("") + ylab("") + 
  scale_color_manual(values = rep("#2171b5", 13)) + 
  scale_x_continuous(labels = lab) + scale_y_continuous(limits = c(0, 1)) + 
  theme_classic() + theme(legend.position = "none") + ggtitle("Fish-Anemone") +
  facet_grid(rows = vars(m), cols = vars(effmeasure)); fane

antp <- spdata %>% filter(Code == "Ant-P") %>% # 168 = 4 * 6 * 7ID
  ggplot(aes(x = Partition, y = value, color = ID)) + 
  geom_line(alpha = 0.2) + geom_point(alpha = 0.2) + xlab("") + ylab("") + 
  scale_color_manual(values = rep("#2171b5", 7)) + 
  scale_x_continuous(labels = lab) + scale_y_continuous(limits = c(0, 1)) + 
  theme_classic() + theme(legend.position = "none") + ggtitle("Ant-Plant") +
  facet_grid(rows = vars(m), cols = vars(effmeasure)); antp

ylab <- textGrob("Proportion of effects", gp = gpar(fontsize = 14), rot = 90)
xlab <- textGrob("Partition", gp = gpar(fontsize = 14))

pa <- plot_grid(paho, phb, herp, prpr, nrow = 2, scale = 1.01)
pa <- grid.arrange(arrangeGrob(pa, left = ylab, bottom = xlab))
#ggsave("../../Dropbox/Kate_Manuscripts/Hierarchical Structure/EL_Submission2/FigSI_SgMSP_Ant.svg", pa, width = 28, height = 20, units = "cm", bg = "white")

pm <- plot_grid(polp, sdp, fane, antp, nrow = 2, scale = 1.01)
pm <- grid.arrange(arrangeGrob(pm, left = ylab, bottom = xlab))
#ggsave("../../Dropbox/Kate_Manuscripts/Hierarchical Structure/EL_Submission2/FigSI_SgMSP_Mut.svg", pm, width = 28, height = 20, units = "cm", bg = "white")

# -------------------------------------------------------------------------
# ---- Fig X: Sector-module partition effect as significance profiles -----

spdata <- partition_df %>% select(ID, Code, effmeasure, m, contains("S")) %>%
  mutate(tot = inSinM + inSbetM + betSinM + betSbetM) %>% 
  mutate(inSinM = inSinM/tot, inSbetM = inSbetM/tot, betSinM = betSinM/tot, betSbetM = betSbetM/tot) %>% 
  select(-tot) %>% mutate(effmeasure = recode(effmeasure, deff = "Direct", ieff = "Indirect")) %>%
  pivot_longer(cols = 5:8, names_to = "Partition", values_to = "value") %>% 
  mutate(Partition = factor(Partition, levels = c("inSinM", "inSbetM", "betSinM", "betSbetM"))) %>%
  mutate(Partition = recode_factor(Partition, "inSinM" = "1", "inSbetM" = "2", "betSinM" = "3", "betSbetM" = "4")) %>% # trick: numeric to plot geom_line
  mutate(Partition = as.numeric(Partition))

lab <- c( "inM\ninSc", "betM\ninSc", "inM\nbetSc", "betM\nbetSc")

paho <- spdata %>% filter(Code == "Pa-Ho") %>% # 1944 = 4 * 6 * 81ID
  ggplot(aes(x = Partition, y = value, color = ID)) + 
  geom_line(alpha = 0.2) + geom_point(alpha = 0.2) + xlab("") + ylab("") + 
  scale_color_manual(values = rep("#d94801", 81)) + 
  scale_x_continuous(labels = lab) + scale_y_continuous(limits = c(0, 1)) + 
  theme_classic() + theme(legend.position = "none") + ggtitle("Parasite-Host") + 
  facet_grid(rows = vars(m), cols = vars(effmeasure)); paho

phb <- spdata %>% filter(Code == "Ph-B") %>% # 912 = 4 * 6 * 38ID
  ggplot(aes(x = Partition, y = value, color = ID)) + 
  geom_line(alpha = 0.2) + geom_point(alpha = 0.2) + xlab("") + ylab("") + 
  scale_color_manual(values = rep("#d94801", 38)) + 
  scale_x_continuous(labels = lab) + scale_y_continuous(limits = c(0, 1)) + 
  theme_classic() + theme(legend.position = "none") + ggtitle("Phage-Bacteria") +
  facet_grid(rows = vars(m), cols = vars(effmeasure)); phb

herp <- spdata %>% filter(Code == "Her-P") %>% # 768 = 4 * 6 * 32ID
  ggplot(aes(x = Partition, y = value, color = ID)) + 
  geom_line(alpha = 0.2) + geom_point(alpha = 0.2) + xlab("") + ylab("") + 
  scale_color_manual(values = rep("#d94801", 32)) + 
  scale_x_continuous(labels = lab) + scale_y_continuous(limits = c(0, 1)) + 
  theme_classic() + theme(legend.position = "none") + ggtitle("Herbivore-Plant") +
  facet_grid(rows = vars(m), cols = vars(effmeasure)); herp

prpr <- spdata %>% filter(Code == "Pr-Pr") %>% # 72 = 4 * 6 * 3ID
  ggplot(aes(x = Partition, y = value, color = ID)) + 
  geom_line(alpha = 0.2) + geom_point(alpha = 0.2) + xlab("") + ylab("") + 
  scale_color_manual(values = rep("#d94801", 3)) + 
  scale_x_continuous(labels = lab) + scale_y_continuous(limits = c(0, 1)) + 
  theme_classic() + theme(legend.position = "none") + ggtitle("Predator-Prey") +
  facet_grid(rows = vars(m), cols = vars(effmeasure)); prpr

polp <- spdata %>% filter(Code == "Pol-P") %>% # 3672 = 4 * 6 * 153ID
  ggplot(aes(x = Partition, y = value, color = ID)) + 
  geom_line(alpha = 0.2) + geom_point(alpha = 0.2) + xlab("") + ylab("") + 
  scale_color_manual(values = rep("#2171b5", 153)) + 
  scale_x_continuous(labels = lab) + scale_y_continuous(limits = c(0, 1)) + 
  theme_classic() + theme(legend.position = "none") + ggtitle("Pollinator-Plant") +
  facet_grid(rows = vars(m), cols = vars(effmeasure)); polp

sdp <- spdata %>% filter(Code == "SD-P") %>% # 936 = 4 * 6 * 39ID
  ggplot(aes(x = Partition, y = value, color = ID)) + 
  geom_line(alpha = 0.2) + geom_point(alpha = 0.2) + xlab("") + ylab("") + 
  scale_color_manual(values = rep("#2171b5", 39)) + 
  scale_x_continuous(labels = lab) + scale_y_continuous(limits = c(0, 1)) + 
  theme_classic() + theme(legend.position = "none") + ggtitle("Seed Disperser-Plant") +
  facet_grid(rows = vars(m), cols = vars(effmeasure)); sdp

fane <- spdata %>% filter(Code == "F-Ane") %>% # 312 = 4 * 6 * 13ID
  ggplot(aes(x = Partition, y = value, color = ID)) + 
  geom_line(alpha = 0.2) + geom_point(alpha = 0.2) + xlab("") + ylab("") + 
  scale_color_manual(values = rep("#2171b5", 13)) + 
  scale_x_continuous(labels = lab) + scale_y_continuous(limits = c(0, 1)) + 
  theme_classic() + theme(legend.position = "none") + ggtitle("Fish-Anemone") +
  facet_grid(rows = vars(m), cols = vars(effmeasure)); fane

antp <- spdata %>% filter(Code == "Ant-P") %>% # 168 = 4 * 6 * 7ID
  ggplot(aes(x = Partition, y = value, color = ID)) + 
  geom_line(alpha = 0.2) + geom_point(alpha = 0.2) + xlab("") + ylab("") + 
  scale_color_manual(values = rep("#2171b5", 7)) + 
  scale_x_continuous(labels = lab) + scale_y_continuous(limits = c(0, 1)) + 
  theme_classic() + theme(legend.position = "none") + ggtitle("Ant-Plant") +
  facet_grid(rows = vars(m), cols = vars(effmeasure)); antp

pa <- plot_grid(paho, phb, herp, prpr, nrow = 2, scale = 1.01)
pa <- grid.arrange(arrangeGrob(pa, left = ylab, bottom = xlab))
#ggsave("../../Dropbox/Kate_Manuscripts/Hierarchical Structure/EL_Submission2/FigSI_MScSP_Ant.svg", pa, width = 28, height = 20, units = "cm", bg = "white")

pm <- plot_grid(polp, sdp, fane, antp, nrow = 2, scale = 1.01)
pm <- grid.arrange(arrangeGrob(pm, left = ylab, bottom = xlab))
#ggsave("../../Dropbox/Kate_Manuscripts/Hierarchical Structure/EL_Submission2/FigSI_MScSP_Mut.svg", pm, width = 28, height = 20, units = "cm", bg = "white")

# -------------------------------------------------------------------------
# ------------ Fig S8: Plotting effect-link ratio inside groups -----------

points <- partition_df %>% # Summary of eff/int plot
  mutate(m = factor(m, levels = c("0.9", "0.5", "0.1"))) %>% 
  filter(effmeasure == "deff" & m == 0.1 | effmeasure == "ieff") %>% # deff only m = 0.1; all ieff 
  mutate(inS = inSinM + inSbetM, inM = inSinM + betSinM, inMt = inMdinMt + betMdinMt) %>% # raw eff in groups
  mutate(totSM = inSinM + inSbetM + betSinM + betSbetM, # tot raw eff
         totMM = inMdinMt + inMdbetMt + betMdinMt + betMdbetMt) %>%
  mutate(inS = inS/totSM, inM = inM/totSM, inMt = inMt/totMM) %>% 
  select(ID, IntType, effmeasure, m, inS, inM, inMt) %>% 
  left_join(net_struct, by = c("ID")) %>% # add prop of interact in groups
  select(ID, IntType = IntType.y, Code, effmeasure, m, inS, inM, inMt, LinkSct, LinkMod, LinkMot) %>%
  filter(!ID %in% c("Herb_Pearse&Alterm", "M_PL_062")) %>% # ID with NA
  mutate(Sct = inS/LinkSct, Mod = inM/LinkMod, Mot = inMt/LinkMot) %>% # effs per interact
  select(ID, IntType, Code, effmeasure, m, Sct, Mod, Mot) %>% 
  group_by(IntType, effmeasure, m) %>% 
  summarise(Sct = mean(Sct), Mod = mean(Mod), Mot = mean(Mot, na.rm = TRUE)) %>% # NA due to 0
  pivot_longer(4:6, names_to = "group", values_to = "value") %>% 
  mutate(group = factor(group, levels = c("Mot", "Mod", "Sct"))) %>% 
  mutate(group = recode_factor(group, Mot = "Subg"))
lines <- lines_from_points_grp(points)

ggplot(points, aes(x = group, y = value)) + geom_point(aes(color = effmeasure, alpha = m)) +
  geom_segment(data = lines, aes(x = x, y = y, xend = xend, yend = yend, color = effmeasure, alpha = m)) +
  geom_hline(yintercept = 1, size = .5, colour = "grey30", alpha = 0.3) +
  xlab("") + ylab("Effect-link ratio inside groups") +
  facet_wrap(~IntType, nrow = 2) + theme_classic() +
  theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 14), axis.title.x = element_text(size = 14), strip.text = element_text(size = 12)) + scale_color_manual(values = c("#00BFC4", "#F8766D"))

# -------------------------------------------------------------------------
# ------ Fig S9: Plotting relative strength of direct x indirect eff ------

partition_df %>% mutate(sum = inSinM + inSbetM + betSinM + betSbetM) %>% 
  select(ID, IntType, Code, effmeasure, m, sum) %>% 
  pivot_wider(names_from = effmeasure, values_from = sum) %>%
  mutate(Direct = deff/(deff + ieff), Indirect = ieff/(deff + ieff)) %>% # proportion 
  group_by(Code, m) %>% summarise(Direct = mean(Direct), Indirect = mean(Indirect)) %>% # 8Code*3m
  ungroup() %>% pivot_longer(3:4, names_to = "effmeasure", values_to = "value") %>%
  mutate(effmeasure = factor(effmeasure, levels = c("Indirect", "Direct"))) %>% 
  ggplot(aes(x = Code, y = value, fill = effmeasure)) + geom_bar(stat = "identity") +
  labs(x = "", y = "Relative strength of effects") + theme_classic() +
  theme(legend.position = "top", legend.title = element_blank()) + facet_wrap(~m)
