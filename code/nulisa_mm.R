# LIBRARIES -------------------------------------------------------------------------------------------------------
library(dplyr)
library(rio)
library(stringr)
library(lme4)
library(lmeresampler)
library(lmerTest)
library(ggResidpanel)
library(emmeans)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(tidyverse)
library(corrplot)
library(boot)

# DATA FORMATTING -------------------------------------------------------------------------------------------------
# Convert plate to long format
plate1 <- rio::import("/Volumes/Boyle-Lab/Shared_Data/PROJECT_Rux/NUElisa/07-09-24 GHVAP-SI53-Malaria CHIM Study_Plate 1_Report.xlsx") %>% 
  pivot_longer(., cols = c(2:ncol(.)), names_to = "R_Day", values_to = "Response") %>% 
  mutate(Plate = "Plate_1")
plate2 <- rio::import("/Volumes/Boyle-Lab/Shared_Data/PROJECT_Rux/NUElisa/07-15-24 GHVAP-SI53-Malaria CHIM Study_Plate 2_Report.xlsx") %>% 
  pivot_longer(., cols = c(2:ncol(.)), names_to = "R_Day", values_to = "Response") %>% 
  mutate(Plate = "Plate_2")
plate3 <- rio::import("/Volumes/boyle-lab/Shared_Data/PROJECT_Rux/NUElisa/05-21-25 NULISAseq_Inflammation_MichelleBoyle.xlsx") %>% 
  pivot_longer(., cols = c(2:ncol(.)), names_to = "R_Day", values_to = "Response") %>% 
  mutate(Plate = "Plate_3")

# Let's combine
nulisa <- plate1 %>% 
  rbind(plate2) %>% 
  rbind(plate3) %>% 
  separate(., R_Day, c("ID", "Days"), remove = FALSE) %>% 
  filter(., !str_detect(ID, "NC|IPC|CONS|SC")) %>%  # Remove controls
  mutate(Day = as.integer(str_extract(Days, "\\d+"))) %>% 
  mutate(Phase = ifelse(Day < 91, "ACUTE", "REINFECTION")) %>% 
  mutate(Timepoint = case_when(Day == 0 ~ "I",
                               Day %in% c(8,9) ~ "T",
                               Day %in% c(11,12) ~ "T+3",
                               Day %in% c(15,16) ~ "T+7",
                               Day %in% c(28,29) ~ "T+21",
                               Day == 91 ~ "rI",
                               Day == 98 ~ "Pre_rT",
                               Day %in% c(99,100) ~ "rT",
                               ID %in% c("R015","R017") & Day == 101 ~ "rT+3",
                               ID %in% c("R001","R006","R007","R018","R019") & Day == 101 ~ "rT",
                               Day %in% c(102,103,104) ~ "rT+3",
                               Day %in% c(106,107,108) ~ "rT+7",
                               Day > 110 ~ "rT+21",
                               TRUE ~ NA))
  # pivot_wider(., names_from = "targetName", values_from = "Response")

##read in the D/R and treatment groups
RandD_treatment <- read.csv("/Volumes/Boyle-Lab/Shared_Data/PROJECT_Rux/NUElisa/RandD_treatment.csv") %>% 
  select(., -c(1)) %>% 
  rename(ID = "R_Number", SampleID = "Sample")

nulisa2 <- left_join(nulisa, RandD_treatment, by="ID")

nulisa2 <- nulisa2 %>% 
  pivot_wider(names_from = targetName, values_from = Response)




# FUNCTION --------------------------------------------------------------------------------------------------------
source("/functions/do.mm_Sep25.R")

# Export xlsx
write.xlsx(nulisa2, "/Users/damian.oyong/Library/CloudStorage/OneDrive-BurnetInstitute/Burnet/Projects/Rux/Manuscript/Clinical/SData/Nulisa.xlsx")

# Run analysis ----------------------------------------------------------------------------------------------------
unique(nulisa2$Timepoint)

# First Infection -------------------------------------------------------------------------------------------------

# I vs T
setwd("/Users/damian.oyong/Library/CloudStorage/OneDrive-BurnetInstitute/Burnet/Projects/Rux/NULISA/Mixed_model/IvT")

IvT <- filter(nulisa2, Timepoint %in% c("I", "T")) %>% 
  mutate(Treatment = factor(Treatment, levels = c("PBO","RUX")),
         Timepoint = factor(Timepoint, levels = c("I","T")))

count_data <- IvT %>% 
  select(., c(10:ncol(.)))

meta_data <- IvT %>% 
  select(., -c(10:ncol(.)))

levels(meta_data$Timepoint)
levels(meta_data$Treatment)

nulisa_IvT <- do.mm2(df_count = mm.impute(count_data), df_meta = meta_data,
                     fixed1 = "Timepoint", random1 = "ID")
saveRDS(nulisa_IvT, "nulisa_IvT.rds")

# T vs TP
setwd("/Users/damian.oyong/Library/CloudStorage/OneDrive-BurnetInstitute/Burnet/Projects/Rux/NULISA/Mixed_model/TvTP")

TvTP <- filter(nulisa2, Timepoint %in% c("T", "T+3", "T+7", "T+21")) %>% 
  mutate(Treatment = factor(Treatment, levels = c("PBO","RUX")),
         Timepoint = factor(Timepoint, levels = c("T", "T+3", "T+7", "T+21")))

count_data <- TvTP %>% 
  select(., c(10:ncol(.)))

meta_data <- TvTP %>% 
  select(., -c(10:ncol(.)))

levels(meta_data$Timepoint)
levels(meta_data$Treatment)

nulisa_TvTP <- do.mm(df_count = mm.impute(count_data), df_meta = meta_data, contrast = TRUE,
                     fixed1 = "Timepoint", fixed2 = "Treatment", random1 = "ID")
saveRDS(nulisa_TvTP, "nulisa_TvTP.rds")

