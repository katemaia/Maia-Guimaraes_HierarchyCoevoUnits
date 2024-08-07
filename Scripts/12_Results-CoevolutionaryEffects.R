# -------------- Manuscript results: coevolutionary effects ---------------

# Author: Kate P Maia
# Checked: 08/2024

# Reads data generated in previous scripts to generate manuscript results on coevolutionary effects across groups and the hierarchical structure of ecological networks. 

# --------------------- Loading library, code and data --------------------

library(tidyverse)

net_struct <- read.table("./Outputs/05_Net-Level_Struct.txt", header = T, sep = "\t")

part_SM <- read.table("./Outputs/09_Partition-Effects-SctMd-Dataframe.txt", header = T, sep = "\t")
part_MM <- read.table("./Outputs/10_Partition-Effects-MdMt-Dataframe.txt", header = T, sep = "\t")

# Combining effects dataframes (partSM and partMM)
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

# MD and M differ because in SM they include Pearse and Altermatt and M_PL_062

# -------------------------------------------------------------------------
# -- Results on (in)direct effects and hierarchical structure (Subg/Mod) --

partition_df %>% # MEAN
  mutate(totMM = inMdinMt + inMdbetMt + betMdinMt + betMdbetMt) %>%
  mutate(inMdinMt = inMdinMt/totMM, inMdbetMt = inMdbetMt/totMM, betMdinMt = betMdinMt/totMM, betMdbetMt = betMdbetMt/totMM) %>%
  select(-totMM, -Code) %>% 
  group_by(effmeasure, m) %>% # MEAN
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
  select(-totMM, -Code) %>% 
  group_by(effmeasure, m) %>% # MEAN
  summarise(inMdinMt = sd(inMdinMt, na.rm = TRUE), inMdbetMt = sd(inMdbetMt, na.rm = TRUE), 
            betMdinMt = sd(betMdinMt, na.rm = TRUE), betMdbetMt = sd(betMdbetMt, na.rm = TRUE)) %>% 
  pivot_wider(names_from = "effmeasure", values_from = 3:6) %>% 
  mutate(inMdinMt_D = inMdinMt_deff*100, inMdbetMt_D = inMdbetMt_deff*100, 
         betMdinMt_D = betMdinMt_deff*100, betMdbetMt_D = betMdbetMt_deff*100, 
         inMdinMt_I = inMdinMt_ieff*100, inMdbetMt_I = inMdbetMt_ieff*100,  
         betMdinMt_I = betMdinMt_ieff*100, betMdbetMt_I = betMdbetMt_ieff*100) %>% 
  select(-contains("deff"), -contains("ieff"))

# -------------------------------------------------------------------------
# -- Results on (in)direct effects and hierarchical structure (Table S2) --

mean <- partition_df %>% # Mean
  mutate(totMM = inMdinMt + inMdbetMt + betMdinMt + betMdbetMt) %>%
  mutate(inMdinMt = inMdinMt/totMM, inMdbetMt = inMdbetMt/totMM, betMdinMt = betMdinMt/totMM, betMdbetMt = betMdbetMt/totMM) %>%
  select(-totMM) %>% group_by(Code, effmeasure, m) %>%
  summarise(inMdinMt = mean(inMdinMt, na.rm = TRUE), inMdbetMt = mean(inMdbetMt, na.rm = TRUE), 
            betMdinMt = mean(betMdinMt, na.rm = TRUE), betMdbetMt = mean(betMdbetMt, na.rm = TRUE)) %>% 
  pivot_wider(names_from = "effmeasure", values_from = 4:7) %>% 
  mutate(inMdinMt_D = inMdinMt_deff*100, inMdbetMt_D = inMdbetMt_deff*100, 
         betMdinMt_D = betMdinMt_deff*100, betMdbetMt_D = betMdbetMt_deff*100, 
         inMdinMt_I = inMdinMt_ieff*100, inMdbetMt_I = inMdbetMt_ieff*100,  
         betMdinMt_I = betMdinMt_ieff*100, betMdbetMt_I = betMdbetMt_ieff*100) %>% 
  select(-contains("deff"), -contains("ieff")) %>% ungroup()

sd <- partition_df %>% # Sd
  mutate(totMM = inMdinMt + inMdbetMt + betMdinMt + betMdbetMt) %>%
  mutate(inMdinMt = inMdinMt/totMM, inMdbetMt = inMdbetMt/totMM, betMdinMt = betMdinMt/totMM, betMdbetMt = betMdbetMt/totMM) %>%
  select(-totMM) %>% group_by(Code, effmeasure, m) %>%
  summarise(inMdinMt = sd(inMdinMt, na.rm = TRUE), inMdbetMt = sd(inMdbetMt, na.rm = TRUE), 
            betMdinMt = sd(betMdinMt, na.rm = TRUE), betMdbetMt = sd(betMdbetMt, na.rm = TRUE)) %>% 
  pivot_wider(names_from = "effmeasure", values_from = 4:7) %>% 
  mutate(inMdinMt_D = inMdinMt_deff*100, inMdbetMt_D = inMdbetMt_deff*100, 
         betMdinMt_D = betMdinMt_deff*100, betMdbetMt_D = betMdbetMt_deff*100, 
         inMdinMt_I = inMdinMt_ieff*100, inMdbetMt_I = inMdbetMt_ieff*100,  
         betMdinMt_I = betMdinMt_ieff*100, betMdbetMt_I = betMdbetMt_ieff*100) %>% 
  select(-contains("deff"), -contains("ieff")) %>% ungroup()

mean <- as.data.frame(mean); sd <- as.data.frame(sd)
mean[,-c(1:2)] <- round(mean[,-c(1:2)], digits = 2)
sd[,-c(1:2)] <- round(sd[,-c(1:2)], digits = 2)

data.frame(Code = mean$Code, m = mean$m, # Table S2
           inMdinMt_D = paste(mean$inMdinMt_D, sd$inMdinMt_D, sep = "_"),
           inMdbetMt_D = paste(mean$inMdbetMt_D, sd$inMdbetMt_D, sep = "_"),
           betMdinMt_D = paste(mean$betMdinMt_D, sd$betMdinMt_D, sep = "_"),
           betMdbetMt_D = paste(mean$betMdbetMt_D, sd$betMdbetMt_D, sep = "_"),
           inMdinMt_I = paste(mean$inMdinMt_I, sd$inMdinMt_I, sep = "_"),
           inMdbetMt_I = paste(mean$inMdbetMt_I, sd$inMdbetMt_I, sep = "_"),
           betMdinMt_I = paste(mean$betMdinMt_I, sd$betMdinMt_I, sep = "_"),
           betMdbetMt_I = paste(mean$betMdbetMt_I, sd$betMdbetMt_I, sep = "_"))

# -------------------------------------------------------------------------
# -- Results on (in)direct effects and hierarchical structure (Mod/Sect) --

partition_df %>% # MEAN
  mutate(totSM = inSinM + inSbetM + betSinM + betSbetM) %>%
  mutate(inSinM = inSinM/totSM, inSbetM = inSbetM/totSM, betSinM = betSinM/totSM, betSbetM = betSbetM/totSM) %>%
  select(-totSM, -Code) %>% 
  group_by(effmeasure, m) %>% # MEAN
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
  select(-totSM, -Code) %>% 
  group_by(effmeasure, m) %>% # MEAN
  summarise(inSinM = sd(inSinM), inSbetM = sd(inSbetM), betSinM = sd(betSinM), betSbetM = sd(betSbetM)) %>% 
  pivot_wider(names_from = "effmeasure", values_from = 3:6) %>% 
  mutate(inSinM_D = inSinM_deff*100, inSbetM_D = inSbetM_deff*100, 
         betSinM_D = betSinM_deff*100, betSbetM_D = betSbetM_deff*100, 
         inSinM_I = inSinM_ieff*100, inSbetM_I = inSbetM_ieff*100,  
         betSinM_I = betSinM_ieff*100, betSbetM_I = betSbetM_ieff*100) %>% 
  select(-contains("deff"), -contains("ieff"))

# -------------------------------------------------------------------------
# -- Results on (in)direct effects and hierarchical structure (Table S3) --

mean <- partition_df %>% mutate(totSM = inSinM + inSbetM + betSinM + betSbetM) %>%
  mutate(inSinM = inSinM/totSM, inSbetM = inSbetM/totSM, betSinM = betSinM/totSM, betSbetM = betSbetM/totSM) %>% 
  select(-totSM) %>% group_by(Code, effmeasure, m) %>%
  summarise(inSinM = mean(inSinM), inSbetM = mean(inSbetM), betSinM = mean(betSinM), betSbetM = mean(betSbetM)) %>% 
  pivot_wider(names_from = "effmeasure", values_from = 4:7) %>% 
  mutate(inSinM_D = inSinM_deff*100, inSbetM_D = inSbetM_deff*100, 
         betSinM_D = betSinM_deff*100, betSbetM_D = betSbetM_deff*100, 
         inSinM_I = inSinM_ieff*100, inSbetM_I = inSbetM_ieff*100,  
         betSinM_I = betSinM_ieff*100, betSbetM_I = betSbetM_ieff*100) %>% 
  select(-contains("deff"), -contains("ieff")) %>% ungroup()

sd <- partition_df %>% mutate(totSM = inSinM + inSbetM + betSinM + betSbetM) %>%
  mutate(inSinM = inSinM/totSM, inSbetM = inSbetM/totSM, betSinM = betSinM/totSM, betSbetM = betSbetM/totSM) %>% 
  select(-totSM) %>% group_by(Code, effmeasure, m) %>%
  summarise(inSinM = sd(inSinM), inSbetM = sd(inSbetM), betSinM = sd(betSinM), betSbetM = sd(betSbetM)) %>% 
  pivot_wider(names_from = "effmeasure", values_from = 4:7) %>% 
  mutate(inSinM_D = inSinM_deff*100, inSbetM_D = inSbetM_deff*100, 
         betSinM_D = betSinM_deff*100, betSbetM_D = betSbetM_deff*100, 
         inSinM_I = inSinM_ieff*100, inSbetM_I = inSbetM_ieff*100,  
         betSinM_I = betSinM_ieff*100, betSbetM_I = betSbetM_ieff*100) %>% 
  select(-contains("deff"), -contains("ieff")) %>% ungroup()

mean <- as.data.frame(mean); sd <- as.data.frame(sd)
mean[,-c(1:2)] <- round(mean[,-c(1:2)], digits = 2)
sd[,-c(1:2)] <- round(sd[,-c(1:2)], digits = 2)

data.frame(Code = mean$Code, m = mean$m, # Table S3
           inSinM_D = paste(mean$inSinM_D, sd$inSinM_D, sep = "_"),
           inSbetM_D = paste(mean$inSbetM_D, sd$inSbetM_D, sep = "_"),
           betSinM_D = paste(mean$betSinM_D, sd$betSinM_D, sep = "_"),
           betSbetM_D = paste(mean$betSbetM_D, sd$betSbetM_D, sep = "_"),
           inSinM_I = paste(mean$inSinM_I, sd$inSinM_I, sep = "_"),
           inSbetM_I = paste(mean$inSbetM_I, sd$inSbetM_I, sep = "_"),
           betSinM_I = paste(mean$betSinM_I, sd$betSinM_I, sep = "_"),
           betSbetM_I = paste(mean$betSbetM_I, sd$betSbetM_I, sep = "_"))