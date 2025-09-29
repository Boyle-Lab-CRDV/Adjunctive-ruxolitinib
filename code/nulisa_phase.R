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

##read in the D/R and treatment groups
RandD_treatment <- read.csv("/Volumes/Boyle-Lab/Shared_Data/PROJECT_Rux/NUElisa/RandD_treatment.csv") %>% 
  select(., -c(1)) %>% 
  rename(ID = "R_Number", SampleID = "Sample")

nulisa2 <- left_join(nulisa, RandD_treatment, by="ID")

nulisa2 <- nulisa2 %>% 
  pivot_wider(names_from = targetName, values_from = Response)




# FUNCTION --------------------------------------------------------------------------------------------------------
source("functions/do.phase.R")


# Run analysis ----------------------------------------------------------------------------------------------------
unique(nulisa2$Timepoint)

nulisa_filter <- filter(nulisa2, Timepoint %in% c("I","T","rI","Pre_rT","rT")) %>% 
  mutate(Phase = factor(Phase, levels = c("ACUTE","REINFECTION")),
         Treatment = factor(Treatment, levels = c("PBO","RUX")),
         Timepoint = factor(Timepoint, levels = c("I","T","rI","Pre_rT","rT")))

# This is where Pre_rT is converted to rT
nulisa_filter2 <- nulisa_filter %>% 
  group_by(ID) %>% 
  filter(!(Timepoint == "rT" & any(Timepoint == "Pre_rT"))) %>% 
  mutate(Timepoint = if_else(Timepoint == "Pre_rT", "rT", Timepoint)) %>% 
  ungroup()

# This is where Pre_rT is converted to rT
count_data2 <- nulisa_filter2 %>% 
  select(., c(10:ncol(.))) %>% 
  select(., -c("LTA|LTB","IFNA2")) %>%  # Plate 3 had no these variables
  select(., -c("FGF19","IFNW1"))  # isSingular

meta_data2 <- nulisa_filter2 %>% 
  select(., -c(10:ncol(.)))

levels(meta_data2$Day)
levels(meta_data2$Phase)
levels(meta_data2$Timepoint)
levels(meta_data2$Treatment)

# Just checking what is the minimum non-zero values
min_df <- nulisa_filter %>% 
  select(c(1,10:ncol(.))) %>% 
  pivot_longer(cols = c(2:ncol(.)), names_to = "Response", values_to = "Values") %>% 
  group_by(Response) %>% 
  summarise(min_value = min(Values, na.rm = TRUE),
            min_non_zero = min(Values[Values > 0], na.rm = TRUE),
            mean_value = mean(Values, na.rm = TRUE)
            )

setwd("/Users/damian.oyong/Library/CloudStorage/OneDrive-BurnetInstitute/Burnet/Projects/Rux/NULISA/Phase/Phase_selected")

phase_selected <- do.phase(df_count = mm.impute(count_data2), df_meta = meta_data2)
saveRDS(phase_selected, "phase_sel.rds")

# Re run with model = "lm" for those with random-effects == 0 or isSingular
count_data2 <- nulisa_filter2 %>% 
  select(., c(10:ncol(.))) %>% 
  select(., c("FGF19","IFNW1")) # isSingular
meta_data2 <- nulisa_filter2 %>% 
  select(., -c(10:ncol(.)))

phase_sel_re <- do.phase(df_count = mm.impute(count_data2), df_meta = meta_data2, model = "lm")
saveRDS(phase_sel_re, "phase_sel_re.rds")

# Combine mm and lm results
phase_sel.comb <- phase_selected
phase_sel.comb$mm_res <- rbind(phase_selected$mm_res, phase_sel_re$mm_res) %>% 
  mutate(fdr = p.adjust(pval, method = "fdr"))
phase_sel.comb$mm_pred <- c(phase_selected$mm_pred, phase_sel_re$mm_pred)
setwd("/Users/damian.oyong/Library/CloudStorage/OneDrive-BurnetInstitute/Burnet/Projects/Rux/NULISA/Phase/Phase_selected")
saveRDS(phase_sel.comb, "phase_sel.comb.rds")



