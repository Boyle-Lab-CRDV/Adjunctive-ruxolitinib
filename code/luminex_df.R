# RUX SUPERNATANT DATA FORMATTING
# 10 JUNE 2025
# R001 and R002 samples were run twice by Stanford team to obtain optimal dilutions, we are using the second run Jun22 for analysis

# IMPORT ----------------------------------------------------------------------------------------------------------
library(openxlsx)

folder_path <- "/Volumes/boyle-lab/Shared_Data/PROJECT_Rux/Luminex/RAW_uRBC_pRBC_supernatant/"

# List all .xlsm files
xlsm_files <- list.files(folder_path, pattern = "\\.xlsm$", full.names = TRUE)

# Function to read "RESULT" sheet from one file
read_result_sheet <- function(file) {
  wb <- loadWorkbook(file)
  sheet_name <- grep("^RESULT$", names(wb), value = TRUE)
  if (length(sheet_name) == 0) return(NULL)  # Skip if no "RESULT" sheet
  data <- readWorkbook(wb, sheet = sheet_name)
  # Use 2nd row as column names
  colnames(data) <- as.character(data[1, ])
  data$file <- basename(file)  # add source filename as column
  return(data)
}

# Apply function to all files
results_list <- map(xlsm_files, read_result_sheet)
results_combined <- bind_rows(results_list)

results_format <- results_combined %>% 
  rename(., Well = 1, Name = 2, Type = 3) %>% 
  filter(., !is.na(Well) & !grepl("S.$", Name)) %>%  # ^ = start of string; S = the character "S"; . = any single character (so total 2 characters); $ asserts the end of the string
  mutate(across(c(4:(ncol(.)-1)), ~parse_number(.)))


# Add variables

no_prefix <- c(
  "CHIM Study (CIRHIS)-H76-H80Panel1-H48_SI53 CHIM Study-1.xlsm",
  "CHIM Study (CIRHIS)-H76-H80Panel1-Plate 2-1.xlsm",
  ""
)

supernatant <- results_format %>% 
  mutate(
    R_Number = case_when(
      str_detect(Name, "R\\d{3}") ~ str_extract(Name, "R\\d{3}"),
      TRUE ~ NA_character_
    ),
    SampleID = case_when(
      str_detect(Name, "D\\d{3}") ~ str_extract(Name, "D\\d{3}"),
      TRUE ~ NA_character_
    ),
    Stim = case_when(
      str_detect(Name, "uRBC") ~ "uRBC",
      str_detect(Name, "pRBC") ~ "pRBC",
      TRUE ~ NA_character_
    ),
    Day = case_when(
      str_detect(Name, "(?i)Day[ -]?(\\d+)") ~ str_extract(Name, "(?i)(?<=Day[ -]?)\\d+"),
      str_detect(Name, "(?<=\\b)\\d+(?=[- ]uRBC|[- ]pRBC)") ~ str_extract(Name, "(?<=\\b)\\d+(?=[- ]uRBC|[- ]pRBC)"),
      TRUE ~ NA_character_
    )
  ) %>% 
  filter(., !(is.na(R_Number) & is.na(SampleID)))
  
RandD_treatment <- read.csv("/Volumes/Boyle-Lab/Shared_Data/PROJECT_Rux/NUElisa/RandD_treatment.csv") %>% 
  select(., -c(1)) %>% 
  rename(SampleID = "Sample")

# Create lookup vectors
rnum_by_sampleid <- setNames(RandD_treatment$R_Number, RandD_treatment$SampleID)
sampleid_by_rnum <- setNames(RandD_treatment$SampleID, RandD_treatment$R_Number)

supernatant_full <- supernatant %>%
  mutate(
    # Fill missing R_Number using SampleID if available
    R_Number = if_else(
      is.na(R_Number) & !is.na(SampleID),
      rnum_by_sampleid[SampleID],
      R_Number
    ),
    # Fill missing SampleID using R_Number if available
    SampleID = if_else(
      is.na(SampleID) & !is.na(R_Number),
      sampleid_by_rnum[as.character(R_Number)],
      SampleID
    ),
    Day = as.numeric(Day)
  ) %>% 
  left_join(., RandD_treatment) %>% 
  select(-c(Type, Well))

setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # This only works if you are using RStudio
saveRDS(supernatant_full, "supernatant_formatted_4Sep25.rds")

