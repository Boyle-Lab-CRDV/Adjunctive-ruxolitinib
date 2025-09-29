#---------------------------------------------------- LIBRARY ----------------------------------------------------
library(lme4)
library(ggResidpanel)
library(lmeresampler)
library(dplyr)
library(tidyr)
library(gtools)
library(openxlsx)
library(cowplot)
library(boot)
library(emmeans)
library(tibble)
library(ggplot2)
library(ggpubr)
library(stringr)

#---------------------------------------------------- INPUT ----------------------------------------------------
# Set your dataframe in a wide format, where each Response variable is in its own column
# IMPORTANT: Factorise your fixed effect variables, use factor(df, levels = c())
# With mixed model, all factors within a fixed effect variable will be compared to the first factor.
# E.g. if levels(meta$fixed1) = "A" "B" "C", output will show comparisons vs "A"


setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # This only works if you are using RStudio

### Import datafile
df <- readRDS("/Volumes/boyle-lab/Shared_Data/PROJECT_Rux/Luminex/RAW_uRBC_pRBC_supernatant/supernatant_formatted_4Sep25.rds") %>% 
  select(., c(Name, R_Number, SampleID, Stim, Day, Treatment), everything()) %>% 
  rename_with(~ str_replace_all(., "[/-]", "_")) %>%     # replace / or - with _
  rename_with(~ str_replace_all(., "\\(LTA\\)", "")) %>% # remove (LTA)
  rename_with(~ str_replace_all(., "__+", "_"))          # collapse any double underscores if created

pt_data <- rio::import("/Users/damian.oyong/Library/CloudStorage/OneDrive-BurnetInstitute/Burnet/Projects/Rux/Luminex/patient_data.xlsx") %>% 
  select(ID, Day, Timepoint) %>% 
  rename(SampleID = ID)

prbc <- df %>% 
  filter(., Stim == "pRBC") %>% 
  mutate(IFN.IL10 = IFNG / IL10,
         IFN.IL10 = if_else(is.infinite(IFN.IL10) | is.nan(IFN.IL10), 0, IFN.IL10),
         Timepoint = case_when(Day == 0 ~ "I",
                               Day %in% c(8,9) ~ "T",
                               Day %in% c(11,12) ~ "T+3",
                               Day %in% c(15,16) ~ "T+7",
                               Day %in% c(28,29) ~ "T+21",
                               TRUE ~ NA)) %>% 
  mutate(Treatment = factor(Treatment, levels = c("PBO","RUX")),
         Timepoint = factor(Timepoint, levels = c("I","T","T+3","T+7","T+21")))

prbc2 <- df %>% 
  filter(., Stim == "pRBC" & Day >90) %>% 
  mutate(drel = Day - 91) %>% 
  left_join(pt_data) %>% 
  mutate(Timepoint = case_when(drel %in% c(27,28,29,30) ~ "rT+21",
                               Day == 104 & SampleID != "D030" ~ "rT+3",
                               Day == 104 & SampleID == "D030" ~ "rT+7",
                               Day == 108 ~ "rT+7",
                               TRUE ~ Timepoint)) %>% 
  mutate(IFN.IL10 = IFNG / IL10,
         IFN.IL10 = if_else(is.infinite(IFN.IL10) | is.nan(IFN.IL10), 0, IFN.IL10)) %>% 
  mutate(Treatment = factor(Treatment, levels = c("PBO","RUX")),
         Timepoint = factor(Timepoint, levels = c("rT","rI","rT+3","rT+7","rT+21")))

urbc <- df %>% 
  filter(., Stim == "uRBC") %>% 
  mutate(IFN.IL10 = IFNG / IL10,
         IFN.IL10 = if_else(is.infinite(IFN.IL10) | is.nan(IFN.IL10), 0, IFN.IL10),
         Timepoint = case_when(Day == 0 ~ "I",
                               Day %in% c(8,9) ~ "T",
                               Day %in% c(11,12) ~ "T+3",
                               Day %in% c(15,16) ~ "T+7",
                               Day %in% c(28,29) ~ "T+21",
                               TRUE ~ NA)) %>% 
  mutate(Treatment = factor(Treatment, levels = c("PBO","RUX")),
         Timepoint = factor(Timepoint, levels = c("I","T","T+3","T+7","T+21")))

urbc2 <- df %>% 
  filter(., Stim == "uRBC" & Day >90) %>% 
  mutate(drel = Day - 91) %>% 
  left_join(pt_data %>% select(., c("SampleID","Day","Timepoint"))) %>% 
  mutate(Timepoint = case_when(drel %in% c(27,28,29,30) ~ "rT+21",
                               Day == 104 & SampleID != "D030" ~ "rT+3",
                               Day == 104 & SampleID == "D030" ~ "rT+7",
                               Day == 108 ~ "rT+7",
                               TRUE ~ Timepoint)) %>% 
  mutate(IFN.IL10 = IFNG / IL10,
         IFN.IL10 = if_else(is.infinite(IFN.IL10) | is.nan(IFN.IL10), 0, IFN.IL10)) %>% 
  mutate(Treatment = factor(Treatment, levels = c("PBO","RUX")),
         Timepoint = factor(Timepoint, levels = c("rT","rI","rT+3","rT+7","rT+21")))

df_sub <- df %>%
  filter(., Day < 90) %>% 
  select(-c("Name")) %>% 
  pivot_wider(
    id_cols = c(SampleID, R_Number, Day, Treatment),
    names_from = Stim,
    values_from = where(is.numeric) & !all_of(c("Day"))
  ) %>%
  transmute(
    SampleID, R_Number, Day, Treatment,
    across(
      ends_with("_pRBC"),
      ~ pmax(. - get(sub("_pRBC$", "_uRBC", cur_column())), 0),
      .names = "{sub('_pRBC$', '', .col)}"
    )
  ) %>% 
  mutate(IFN.IL10 = IFNG / IL10,
         IFN.IL10 = if_else(is.infinite(IFN.IL10) | is.nan(IFN.IL10), 0, IFN.IL10),
         Timepoint = case_when(Day == 0 ~ "I",
                               Day %in% c(8,9) ~ "T",
                               Day %in% c(11,12) ~ "T+3",
                               Day %in% c(15,16) ~ "T+7",
                               Day %in% c(28,29) ~ "T+21",
                               TRUE ~ NA)) %>% 
  mutate(Treatment = factor(Treatment, levels = c("PBO","RUX")),
         Timepoint = factor(Timepoint, levels = c("I","T","T+3","T+7","T+21")))


df_sub2 <- df %>%
  filter(., Day > 90) %>% 
  select(-c("Name")) %>% 
  pivot_wider(
    id_cols = c(SampleID, R_Number, Day, Treatment),
    names_from = Stim,
    values_from = where(is.numeric) & !all_of(c("Day"))
  ) %>%
  transmute(
    SampleID, R_Number, Day, Treatment,
    across(
      ends_with("_pRBC"),
      ~ pmax(. - get(sub("_pRBC$", "_uRBC", cur_column())), 0),
      .names = "{sub('_pRBC$', '', .col)}"
    )
  ) %>% 
  mutate(drel = Day - 91) %>% 
  left_join(pt_data %>% select(., c("SampleID","Day","Timepoint"))) %>% 
  mutate(Timepoint = case_when(drel %in% c(27,28,29,30) ~ "rT+21",
                               Day == 104 & SampleID != "D030" ~ "rT+3",
                               Day == 104 & SampleID == "D030" ~ "rT+7",
                               Day == 108 ~ "rT+7",
                               TRUE ~ Timepoint)) %>% 
  mutate(IFN.IL10 = IFNG / IL10,
         IFN.IL10 = if_else(is.infinite(IFN.IL10) | is.nan(IFN.IL10), 0, IFN.IL10)) %>% 
  mutate(Treatment = factor(Treatment, levels = c("PBO","RUX")),
         Timepoint = factor(Timepoint, levels = c("rT","rI","rT+3","rT+7","rT+21"))) 

saveRDS(df_sub, "df_sub.rds")
saveRDS(df_sub2, "df_sub2.rds")

# Export as excel
write.xlsx(df_sub, "/Users/damian.oyong/Library/CloudStorage/OneDrive-BurnetInstitute/Burnet/Projects/Rux/Manuscript/Clinical/SData/Luminex.xlsx")


#---------------------------------------------------- FUNCTION ----------------------------------------------------
source("functions/do.mm_Sep25.R")


##### Main Script #####

# Analyse data I vs T ---------------------------------------------------------------------------------------------
setwd("/Users/damian.oyong/Library/CloudStorage/OneDrive-BurnetInstitute/Burnet/Projects/Rux/Luminex/Supernatant/I_T")

## Split main df into 2 types; count df containing dependent variables and metadata df
df_sub.IT <- df_sub %>% 
  filter(Day < 10) %>% 
  mutate(
         Timepoint = case_when(Day == 0 ~ "I",
                               Day %in% c(8,9) ~ "T",
                               TRUE ~ NA)) %>% 
  mutate(Treatment = factor(Treatment, levels = c("PBO","RUX")),
         Timepoint = factor(Timepoint, levels = c("I","T"))
         )

count_data.IT <- df_sub.IT %>% 
  select(., c(6:ncol(.)-1)) %>%
  select(., -c("IL15","IL18","RANTES_CCL5")) %>%  # IL15,IL18,RANTES_CCL5 are all 0
  select(., -c("IL12P70")) # Need to be analysed with lm.opt = TRUE

meta_data.IT <- df_sub.IT %>% 
  select(., -c(6:ncol(.)-1))

levels(meta_data.IT$Timepoint)
levels(meta_data.IT$Treatment)

super.IT <- do.mm2(df_count = mm.transform(count_data.IT), df_meta = meta_data.IT, lm.opt = FALSE,
                  fixed1 = "Timepoint", random1 = "SampleID")
saveRDS(super.IT, file = "Supernatant_IvT.rds")
super.IT <- readRDS("Supernatant_IvT.rds")

count_data.IT <- df_sub.IT %>% 
  select(., c("IL12P70"))

meta_data.IT <- df_sub.IT %>% 
  select(., -c(6:ncol(.)-1))

levels(meta_data.IT$Timepoint)
levels(meta_data.IT$Treatment)

super.IT.other <- do.mm2(df_count = mm.transform(count_data.IT), df_meta = meta_data.IT, lm.opt = TRUE,
                         fixed1 = "Timepoint", random1 = "SampleID")
saveRDS(super.IT.other, file = "Supernatant.other_IvT.rds")
super.IT.other <- readRDS("Supernatant.other_IvT.rds")

super.IT.comb <- super.IT

super.IT.comb$mm_res <- super.IT$mm_res %>% 
  bind_rows(., super.IT.other$mm_res) %>% 
  dplyr::group_by(term, fixed) %>%
  mutate(fdr = p.adjust(pval, method = "fdr"))
  
super.IT.comb$mm_pred <- c(super.IT$mm_pred, super.IT.other$mm_pred)
saveRDS(super.IT.comb, file = "Supernatant.comb_IvT.rds")


# pRBC
prbc_IT <- prbc %>% 
  filter(Day < 10) %>% 
  mutate(IFN.IL10 = IFNG/IL10,
         Timepoint = case_when(Day == 0 ~ "I",
                               Day %in% c(8,9) ~ "T",
                               TRUE ~ NA)) %>% 
  mutate(Treatment = factor(Treatment, levels = c("PBO","RUX")),
         Timepoint = factor(Timepoint, levels = c("I","T"))
  )

# Analyse data T vs TP ---------------------------------------------------------------------------------------------
setwd("/Users/damian.oyong/Library/CloudStorage/OneDrive-BurnetInstitute/Burnet/Projects/Rux/Luminex/Supernatant/T_TP")

## Split main df into 2 types; count df containing dependent variables and metadata df
df_sub.TP <- df_sub %>% 
  filter(Day < 90 & Day > 0) %>% 
  mutate(
         Timepoint = case_when(Day %in% c(8,9) ~ "T",
                               Day %in% c(11,12) ~ "T+3",
                               Day %in% c(15,16) ~ "T+7",
                               Day %in% c(28,29) ~ "T+21",
                               TRUE ~ NA)) %>% 
  mutate(Treatment = factor(Treatment, levels = c("PBO","RUX")),
         Timepoint = factor(Timepoint, levels = c("T","T+3","T+7","T+21")))

count_data.TP <- df_sub.TP %>% 
  select(., c(6:ncol(.)-1)) %>% 
  select(., -c("IL15","IL18","RANTES_CCL5")) # These are all zeroes

meta_data.TP <- df_sub.TP %>% 
  select(., c(R_Number, SampleID, Day, Treatment, Timepoint))

levels(meta_data.TP$Timepoint)
levels(meta_data.TP$Treatment)

super.TP <- do.mm(df_count = mm.transform(count_data.TP), df_meta = meta_data.TP, contrast = TRUE,
                  intercept = TRUE, fixed1 = "Timepoint", fixed2 = "Treatment", random1 = "SampleID")
saveRDS(super.TP, file = "Supernatant_TvTP.rds")
super.TP <- readRDS("Supernatant_TvTP.rds")




# Analyse data rT vs rTP ------------------------------------------------------------------------------------------
setwd("/Users/damian.oyong/Library/CloudStorage/OneDrive-BurnetInstitute/Burnet/Projects/Rux/Luminex/Supernatant/rTP")

## Split main df into 2 types; count df containing dependent variables and metadata df
count_data.rTP <- df_sub2 %>% 
  select(., c(6:ncol(.))) %>% 
  select(., -c(EGF,IL3,IL15,IL18,RANTES_CCL5,TGFA)) %>% 
  select(., -c(drel, Timepoint))
meta_data.rTP <- df_sub2 %>% 
  select(., c(1:4,54))

levels(meta_data.rTP$Timepoint)
levels(meta_data.rTP$Treatment)

super.rTP <- do.mm(df_count = mm.transform(count_data.rTP), df_meta = meta_data.rTP, contrast = TRUE,
                  fixed1 = "Timepoint", fixed2 = "Treatment", random1 = "SampleID")
saveRDS(super.rTP, file = "Supernatant_rTP.rds")
super.rTP <- readRDS("Supernatant_rTP.rds")




# Line Plot 1st Phase -------------------------------------------------------------------------------------------------------

plot.dir <- "/Users/damian.oyong/Library/CloudStorage/OneDrive-BurnetInstitute/Burnet/Projects/Rux/Luminex/Supernatant/Figures/"
dessert <- c("#5C62D6", "#BF2C34")

to.plot <- c("IFNG","IL10","IFN.IL10")


## BG Subtracted ------------------------------------------------------------------------------------------------------------

for(i in unique(super.TP$mm_res$Response)){
  df.plot <- df_sub %>% 
    select(., -c(5:(ncol(.)-1))) %>%
    cbind(.,
          df_sub %>%
            select(., c(5:(ncol(.)-1))) %>%
            mm.transform(., log = "log10"))
  
  y_pos <- rstatix::get_y_position(df.plot, 
                                   formula =  as.formula(paste0(i, "~ Timepoint")),
                                   step.increase = 0.28, fun = "max")
  low = y_pos$y.position[1]
  mid = y_pos$y.position[2]
  up = y_pos$y.position[3]
  
  
  # p values for I vs T
  p_manual <- super.IT$mm_res %>% 
    select(term, Response, pval, fixed) %>% 
    pivot_wider(names_from = fixed, values_from = pval) %>% 
    dplyr::rename(., group2 = term, .y. = Response, Timepoint = c(3)) %>% 
    mutate(group1 = "I", 
           Timepoint = case_when(Timepoint < 0.001 ~ "<0.001",
                                 Timepoint < 0.10 & Timepoint >= 0.001 ~ as.character(round(Timepoint, digits = 3)),
                                 TRUE ~ as.character(round(Timepoint, digits = 2))),
           group2 = factor(group2, levels = c("I","T"))
    ) %>% 
    group_by(group2) %>% 
    mutate(xmin = as.numeric(group2), xmax = as.numeric(group2)) %>% 
    ungroup()
  
  # p values T vs TP
  ## For TPPBO and TX
  p_manual2 <- super.TP$mm_res %>% 
    select(term, Response, pval, fixed) %>% 
    pivot_wider(names_from = fixed, values_from = pval) %>% 
    dplyr::rename(., group2 = term, .y. = Response, Timepoint=c(3)) %>% 
    mutate(group1 = "T", 
           Timepoint = case_when(Timepoint < 0.001 ~ "<0.001",
                                 Timepoint < 0.10 & Timepoint >= 0.001 ~ as.character(round(Timepoint, digits = 3)),
                                 TRUE ~ as.character(round(Timepoint, digits = 2))),
           Treatment = case_when(Treatment < 0.001 ~ "<0.001",
                                 Treatment < 0.10 & Treatment >= 0.001 ~ as.character(round(Treatment, digits = 3)),
                                 TRUE ~ as.character(round(Treatment, digits = 2))),
           group2 = factor(group2, levels = c("T","T+3", "T+7", "T+21"))
    ) %>% 
    group_by(group2) %>% 
    mutate(xmin = as.numeric(group2)+1, xmax = as.numeric(group2)+1) %>% 
    ungroup()
  
  ## For TX
  p_manual3 <- super.TP$mm_contrast %>% 
    select(term, Response, pvalue) %>% 
    dplyr::rename(., group2 = "term", Timepoint = "pvalue", .y. = "Response") %>%  
    mutate(group1 = "T", 
           Timepoint = case_when(Timepoint < 0.001 ~ "<0.001",
                                 Timepoint < 0.10 & Timepoint >= 0.001 ~ as.character(round(Timepoint, digits = 3)),
                                 TRUE ~ as.character(round(Timepoint, digits = 2))),
           group2 = factor(group2, levels = c("RUX_T+3vT", "RUX_T+7vT", "RUX_T+21vT"))
    ) %>% 
    group_by(group2) %>% 
    mutate(xmin = as.numeric(group2)+2, xmax = as.numeric(group2)+2) %>% 
    ungroup()
  
  print(
    ggplot(df.plot, aes(x = Timepoint, y = .data[[i]], group = SampleID))+
      geom_line(data = filter(df.plot, Timepoint %in% c("I","T")), color = "grey", alpha=0.4)+
      geom_line(data = filter(df.plot, !Timepoint %in% c("I")), aes(color = Treatment), alpha=0.4)+
      ## This is to generate the ribbon of SE
      ### For I vs T
      stat_summary(data = super.IT$mm_pred[[i]],
                   aes(x = Timepoint, y = predict),
                   fun.data = function(y) {
                     mean_y <- mean(y)
                     var_predict <- var(super.IT$mm_pred[[i]]$predict)
                     se_y <- sqrt(mean(super.IT$mm_pred[[i]]$se^2) + var_predict)
                     return(c(y = mean_y, ymin = mean_y - se_y, ymax = mean_y + se_y))
                   },
                   geom = "ribbon", fill = "grey", alpha = 0.03) +
      ### For T vs TP
      stat_summary(data = super.TP$mm_pred[[i]],
                   aes(x = Timepoint, y = predict, fill = Treatment),
                   fun.data = function(y) {
                     mean_y <- mean(y)
                     var_predict <- var(super.TP$mm_pred[[i]]$predict)
                     se_y <- sqrt(mean((super.TP$mm_pred[[i]]$se)^2) + var_predict)
                     return(c(y = mean_y, ymin = mean_y - se_y, ymax = mean_y + se_y))
                   },
                   geom = "ribbon", alpha = 0.03) +
      
      ## This is to generate "mean" of predicted values
      ### For I vs T
      stat_summary(data = super.IT$mm_pred[[i]], aes(x = Timepoint, y = predict, colour = "grey", group= "grey"),
                   color = "#4F4F4F",
                   fun = mean, geom = "line", linewidth = 1.5)+
      ### For T vs TP
      stat_summary(data = super.TP$mm_pred[[i]], aes(x = Timepoint, y = predict, color = Treatment, group = Treatment),
                   fun = mean, geom = "line", linewidth = 1.5)+
      
      ## This is for p-values
      ### For I vs T
      stat_pvalue_manual(subset(p_manual, .y.== i), label = "p= {Timepoint}", 
                         y.position = rstatix::get_y_position(filter(df.plot, Timepoint %in% c("I","T")), 
                                                              formula =  as.formula(paste0(i, "~ Timepoint")),
                                                              step.increase = 0.28)$y.position, 
                         xmin="group2", xmax="group1")+
      
      ### For T vs TP
      #### TX
      stat_pvalue_manual(subset(p_manual2, .y.== i), label = "{Treatment}",
                         y.position = up, step.increase = 0,
                         tip.length = 0)+
      annotate("segment", y = up,
               x = min(p_manual2$xmin), xend = max(p_manual2$xmax),
               colour = "red3", linewidth = 0.1) +
      annotate("text", y = up, x = min(p_manual2$xmin-0.4),
               label = "TX")+
      
      #### TPPBO
      stat_pvalue_manual(subset(p_manual2, .y.== i), label = "{Timepoint}",
                         y.position = mid, step.increase = 0,
                         tip.length = 0)+
      annotate("segment", y = mid,
               x = min(p_manual2$xmin), xend = max(p_manual2$xmax),
               colour = "blue3", linewidth = 0.1) +
      annotate("text", y = mid, x = min(p_manual2$xmin-0.4),
               label = "PBO")+
      
      #### RUX
      stat_pvalue_manual(subset(p_manual3, .y.== i), label = "{Timepoint}",
                         y.position = low, step.increase = 0,
                         tip.length = 0)+
      annotate("segment", y = low,
               x = min(p_manual3$xmin), xend = max(p_manual3$xmax),
               colour = "blue3", linewidth = 0.1) +
      annotate("text", y = low, x = min(p_manual3$xmin-0.4),
               label = "RUX")+
      
      theme_bw()+
      scale_color_manual(values = dessert)+
      scale_fill_manual(values = dessert) 
  )
  ggsave(filename = paste0(plot.dir,"/",i,".pdf"), width = 6, height = 6)
}


# Line Plot 2nd Phase -------------------------------------------------------------------------------------------------------

plot.dir <- "/Users/damian.oyong/Library/CloudStorage/OneDrive-BurnetInstitute/Burnet/Projects/Rux/Luminex/Supernatant/Figures/"
dessert <- c("#5C62D6", "#BF2C34")

to.plot <- c("IFNG","TNFA","IL10","IFN.IL10")


## BG Subtracted ------------------------------------------------------------------------------------------------------------

for(i in unique(super.rTP$mm_res$Response)){
  df.plot <- df_sub2 %>% 
    select(., -c(6:(ncol(.)-2), ncol(.))) %>%
    cbind(.,
          df_sub2 %>%
            select(., c(6:(ncol(.)-2), ncol(.))) %>%
            mm.transform(., log = "log10")) %>% 
    mutate(., Timepoint = factor(Timepoint, levels = c("rI","rT","rT+3", "rT+7", "rT+21")))
  
  
  y_pos <- rstatix::get_y_position(df.plot, 
                                   formula =  as.formula(paste0(i, "~ Timepoint")),
                                   step.increase = 0.28, fun = "max")
  low = y_pos$y.position[1]
  mid = y_pos$y.position[2]
  up = y_pos$y.position[3]
  
  # p values T vs TP
  ## For TPPBO and TX
  p_manual2 <- super.rTP$mm_res %>% 
    select(term, Response, fdr, fixed) %>% 
    pivot_wider(names_from = fixed, values_from = fdr) %>% 
    dplyr::rename(., group2 = term, .y. = Response, Timepoint=c(3)) %>% 
    mutate(group1 = "rT", 
           Timepoint = case_when(Timepoint < 0.001 ~ "<0.001",
                                 Timepoint < 0.10 & Timepoint >= 0.001 ~ as.character(round(Timepoint, digits = 3)),
                                 TRUE ~ as.character(round(Timepoint, digits = 2))),
           Treatment = case_when(Treatment < 0.001 ~ "<0.001",
                                 Treatment < 0.10 & Treatment >= 0.001 ~ as.character(round(Treatment, digits = 3)),
                                 TRUE ~ as.character(round(Treatment, digits = 2))),
           group2 = factor(group2, levels = c("rI","rT","rT+3", "rT+7", "rT+21"))
    ) %>% 
    group_by(group2) %>% 
    mutate(xmin = as.numeric(group2), xmax = as.numeric(group2)) %>% 
    ungroup()
  
  ## For TX
  p_manual3 <- super.rTP$mm_contrast %>% 
    select(term, Response, fdr) %>% 
    dplyr::rename(., group2 = "term", Timepoint = "fdr", .y. = "Response") %>%  
    mutate(group1 = "rT", 
           Timepoint = case_when(Timepoint < 0.001 ~ "<0.001",
                                 Timepoint < 0.10 & Timepoint >= 0.001 ~ as.character(round(Timepoint, digits = 3)),
                                 TRUE ~ as.character(round(Timepoint, digits = 2))),
           group2 = factor(group2, levels = c("RUX_rIvrT","rT","RUX_rT+3vrT", "RUX_rT+7vrT", "RUX_rT+21vrT"))
    ) %>% 
    group_by(group2) %>% 
    mutate(xmin = as.numeric(group2), xmax = as.numeric(group2)) %>% 
    ungroup()
  
  print(
    ggplot(df.plot, aes(x = Timepoint, y = .data[[i]], group = SampleID))+
      geom_line(data = df.plot, aes(color = Treatment), alpha=0.4)+
      ## This is to generate the ribbon of SE
      ### For T vs TP
      stat_summary(data = super.rTP$mm_pred[[i]],
                   aes(x = Timepoint, y = predict, fill = Treatment),
                   fun.data = function(y) {
                     mean_y <- mean(y)
                     var_predict <- var(super.rTP$mm_pred[[i]]$predict)
                     se_y <- sqrt(mean((super.rTP$mm_pred[[i]]$se)^2) + var_predict)
                     return(c(y = mean_y, ymin = mean_y - se_y, ymax = mean_y + se_y))
                   },
                   geom = "ribbon", alpha = 0.03) +
      
      ## This is to generate "mean" of predicted values
      ### For T vs TP
      stat_summary(data = super.rTP$mm_pred[[i]], aes(x = Timepoint, y = predict, color = Treatment, group = Treatment),
                   fun = mean, geom = "line", linewidth = 1.5)+
      
      
      ### For T vs TP
      #### TX
      stat_pvalue_manual(subset(p_manual2, .y.== i), label = "{Treatment}",
                         y.position = up, step.increase = 0,
                         tip.length = 0)+
      annotate("segment", y = up,
               x = min(p_manual2$xmin), xend = max(p_manual2$xmax),
               colour = "red3", linewidth = 0.1) +
      annotate("text", y = up, x = min(p_manual2$xmin-0.4),
               label = "TX")+
      
      #### TPPBO
      stat_pvalue_manual(subset(p_manual2, .y.== i), label = "{Timepoint}",
                         y.position = mid, step.increase = 0,
                         tip.length = 0)+
      annotate("segment", y = mid,
               x = min(p_manual2$xmin), xend = max(p_manual2$xmax),
               colour = "blue3", linewidth = 0.1) +
      annotate("text", y = mid, x = min(p_manual2$xmin-0.4),
               label = "PBO")+
      
      #### RUX
      stat_pvalue_manual(subset(p_manual3, .y.== i), label = "{Timepoint}",
                         y.position = low, step.increase = 0,
                         tip.length = 0)+
      annotate("segment", y = low,
               x = min(p_manual3$xmin), xend = max(p_manual3$xmax),
               colour = "blue3", linewidth = 0.1) +
      annotate("text", y = low, x = min(p_manual3$xmin-0.4),
               label = "RUX")+
      
      theme_bw()+
      scale_color_manual(values = dessert)+
      scale_fill_manual(values = dessert) 
  )
  ggsave(filename = paste0(plot.dir,"/REINFECTION_",i,".pdf"), width = 6, height = 6)
}



# Volcano plot ----------------------------------------------------------------------------------------------------
library(ggrepel)

plot.dir

# 1st infection Placebo
super.TP$mm_res %>%
  filter(term != "T" & fixed == "Timepoint" & Response != "IFN.IL10") %>%
  mutate(term = factor(term, levels = c("T+3","T+7","T+21")),
         sig_status = case_when(
           fdr >= 0.05 ~ "nonsig",
           rep.mean > 0 ~ "positive",
           rep.mean < 0 ~ "negative"
         ),
         shape=if_else(pval > 0.05, "ns","signif")
         ) %>%
  ggplot(aes(x = rep.mean, y = -log10(fdr))) +
  geom_point(aes(colour = sig_status, shape = shape), size = 2) +
  geom_text_repel(aes(label = Response)) +
  theme_bw() +
  scale_colour_manual(values = c(
    "positive" = "#EEB422",
    "negative" = "green4",
    "nonsig"   = "grey"
  )) +
  scale_shape_manual(values = c(
    "ns" = 15,
    "signif" = 16
  ))+
  ggtitle("PBO Timepoint") +
  labs(x = "Coefficient", y = "-log10 P-value") +
  facet_wrap(~term, scales = "free") +
  theme(legend.position = "none")
ggsave(filename = paste0(plot.dir,"/Volcano_PBO_TP_fdr.pdf"), width = 12, height = 6)


# 1st infection Rux
super.TP$mm_contrast %>% 
  filter(., Response != "IFN.IL10") %>% 
  mutate(term = factor(term, levels = c("RUX_T+3vT","RUX_T+7vT","RUX_T+21vT")),
         sig_status = case_when(
           fdr >= 0.05 ~ "nonsig",
           rep.mean > 0 ~ "positive",
           rep.mean < 0 ~ "negative"
         ),
         shape=if_else(pvalue > 0.05, "ns","signif")
         ) %>%
  ggplot(., aes(x = rep.mean, y = -log10(fdr)))+
  geom_text_repel(aes(label = Response))+
  geom_point(aes(colour = sig_status, shape = shape), size = 2) +
  theme_bw()+
  scale_colour_manual(values = c(
    "positive" = "#EEB422",
    "negative" = "green4",
    "nonsig"   = "grey"
  )) +
  scale_shape_manual(values = c(
    "ns" = 15,
    "signif" = 16
  ))+
  ggtitle("Rux Timepoint")+
  labs(x = "Coefficient", y = "-log10 P-value")+
  theme(legend.position = "none")+
  facet_wrap(~term, scales = "free")
ggsave(filename = paste0(plot.dir,"/Volcano_RUX_TP_fdr.pdf"), width = 12, height = 6)


# 1st infection TX
super.TP$mm_res %>% 
  filter(., term != "T" & fixed == "Treatment" & Response != "IFN.IL10") %>% 
  mutate(term = factor(term, levels = c("T+3","T+7","T+21")),
         sig_status = case_when(
           fdr >= 0.05 ~ "nonsig",
           rep.mean > 0 ~ "positive",
           rep.mean < 0 ~ "negative"
         ),
         shape=if_else(pval > 0.05, "ns","signif")
         ) %>%
  ggplot(., aes(x = rep.mean, y = -log10(fdr)))+
  geom_text_repel(aes(label = Response))+
  geom_point(aes(colour = sig_status, shape = shape), size = 2) +
  theme_bw()+
  scale_colour_manual(values = c(
    "positive" = "#EEB422",
    "negative" = "green4",
    "nonsig"   = "grey"
  )) +
  scale_shape_manual(values = c(
    "ns" = 15,
    "signif" = 16
  ))+
  ggtitle("Treatment interaction")+
  labs(x = "Coefficient", y = "-log10 P-value")+
  theme(legend.position = "none")+
  facet_wrap(~term, scales = "free")
ggsave(filename = paste0(plot.dir,"/Volcano_TX_fdr.pdf"), width = 12, height = 6)

