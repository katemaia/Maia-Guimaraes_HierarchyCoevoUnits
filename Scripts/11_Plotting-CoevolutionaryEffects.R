# -------------- Manuscript figures: coevolutionary effects ---------------

# Author: Kate P Maia
# Checked: 12/2023

# Reads data generated in previous scripts to plot manuscript figures on coevolutionary effects across groups and the hierarchical structure of ecological networks. 

# --------------------- Loading library, code and data --------------------

library(tidyverse)
library(reshape2)
library(RColorBrewer)
library(GGally)
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

partition_df %>% select(ID, IntType, Code, effmeasure, m, sum) %>% 
  pivot_wider(names_from = effmeasure, values_from = sum) %>%
  mutate(Direct = deff/(deff + ieff), Indirect = ieff/(deff + ieff)) %>% # proportion 
  group_by(Code, m) %>% summarise(Direct = mean(Direct), Indirect = mean(Indirect)) %>% # 8Code*3m
  ungroup() %>% pivot_longer(3:4, names_to = "effmeasure", values_to = "value") %>%
  mutate(effmeasure = factor(effmeasure, levels = c("Indirect", "Direct"))) %>% 
  ggplot(aes(x = Code, y = value, fill = effmeasure)) + geom_bar(stat = "identity") +
  labs(x = "", y = "Relative strength of effects") + theme_classic() +
  theme(legend.position = "top", legend.title = element_blank()) + facet_wrap(~m)