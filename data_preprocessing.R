library(tidyverse)
library(Hmisc)
library(corrplot)
library(mice)
library(naniar)
library(gtsummary)

data_clin <- read.csv("Copia di clinical_variables_MM.csv", sep = ";", na = "#N/D")
data_gen  <- read.csv("CNA_MM.csv", sep = ";", na = "#N/D")

data_gen <- data_gen %>%
  mutate(across(where(is.character), ~tolower(.))) %>%
  mutate(across(-c(1, 2), as.factor))

data_clin <- data_clin %>%
  mutate(across(where(is.character), ~tolower(.))) %>%
  mutate(across(where(is.character), ~na_if(.x, "nv")))

data_clin <- data_clin %>%
  mutate(
    SEX = as.factor(SEX),
    FISH_T_11_14 = as.factor(FISH_T_11_14),
    CREATININE = as.numeric(str_replace(CREATININE, ",", ".")),
    Clearence = as.numeric(str_replace(Clearence, ",", ".")),
    Clearance_m_50 = as.factor(Clearance_m_50),
    HB = as.numeric(str_replace(HB, ",", ".")),
    HB_m_105 = as.factor(HB_m_105),
    PLT = as.numeric(PLT),
    PLT_m_150 = as.factor(PLT_m_150),
    PC_TOT = as.numeric(PC_TOT),
    PC_M_60 = as.factor(PC_M_60),
    R_ISS = as.factor(R_ISS),
    IG_ISOTYPE = as.factor(IG_ISOTYPE),
    LIGHT_CHAIN = as.factor(LIGHT_CHAIN),
    ECOG = as.factor(ECOG),
    CALCIUM = as.numeric(str_replace(CALCIUM, ",", ".")),
    Calcio_M_105 = as.factor(Calcio_M_105),
    PCR = as.numeric(str_replace(PCR, ",", ".")),
    PCR_m_05 = as.factor(PCR_m_05),
    
    PET_POS_NEG = as.factor(PET_POS_NEG),
    BM_DS = as.factor(BM_DS),
    BM_DS_M_4 = as.factor(BM_DS_M_4),
    BM_Suv_max = as.numeric(BM_Suv_max),
    FL_sec_Impetus = as.factor(FL_sec_Impetus),
    FL_DS = as.factor(FL_DS),
    FL_DS_M_4 = as.factor(FL_DS_M_4),
    FL_SUV_max = as.numeric(FL_SUV_max),
    PM = as.factor(PM),
    PM_DS = as.factor(PM_DS),
    PM_SUVmax = as.numeric(PM_SUVmax),
    EM = as.factor(EM),
    EM_DS = as.factor(EM_DS),
    EM_SUV_max = as.numeric(EM_SUV_max),
    
  )

vars_clinical <- c("MPC", "AGE", "SEX", "FISH_T_11_14", "CREATININE", "Clearence", 
                   "HB", "PLT", "PC_TOT", "R_ISS", "IG_ISOTYPE", "LIGHT_CHAIN", 
                   "ECOG", "CALCIUM", "PCR")

vars_imaging <- c("gDNA_tumor_fraction", "cfDNA_tumor_fraction", 
                  "BM_Suv_max", "FL_SUV_max", 
                 "PM_SUVmax", "EM_SUV_max")

vars_survival <- c("PFS_I_EVENT", "PFS_I_MONTHS")

data_clin <- data_clin %>%
  select(all_of(c(vars_imaging, vars_clinical, vars_survival)))


data_clin_clean <- data_clin %>%
  #Exclude patients with no outcome
  drop_na(PFS_I_MONTHS, PFS_I_EVENT) %>%
  # Exclude variables with >50% missingness
  select(where(~ mean(is.na(.)) <= 0.4))

# check missing data
missing_analysis_pre <- naclus(data_clin_clean %>% select(-MPC))

par(mfrow = c(1, 1))
plot(missing_analysis_pre)
naplot(missing_analysis_pre, which = "all")

# Exclusion of patients with highly incomplete clinical dataset
extreme_na_patients <- data_clin_clean %>% 
  filter(rowSums(is.na(select(., -MPC))) >= 6)
data_clin_final <- data_clin_clean %>%
  filter(rowSums(is.na(select(., -MPC))) < 6)

# Clinical Categorical Variables Frequencies
print("--- CLINICAL CATEGORICAL VARIABLES FREQUENCIES ---")
lapply(data_clin_final %>% select(where(is.factor)), table)

# Exclusion of the single non-secretory ("ns") light chain patient
id_ns_light_chain <- data_clin_final %>%
  filter(LIGHT_CHAIN == "ns") %>%
  pull(MPC)
data_clin_final <- data_clin_final %>%
  filter(is.na(LIGHT_CHAIN) | !LIGHT_CHAIN %in% "ns") %>%
  mutate(LIGHT_CHAIN = droplevels(LIGHT_CHAIN))

table(data_clin_final$LIGHT_CHAIN, useNA = "ifany")

data_clin_final <- data_clin_final %>%
  mutate(
    IG_ISOTYPE = case_when(
      IG_ISOTYPE %in% c("igg", "igg/bj") ~ "IgG",
      IG_ISOTYPE %in% c("iga", "iga/bj", "iga+bj") ~ "IgA",
      IG_ISOTYPE == "bj" ~ "BJ",
      IG_ISOTYPE %in% c("igd", "ns") ~ "other",
      TRUE ~ as.character(IG_ISOTYPE)
    ),
    IG_ISOTYPE = droplevels(factor(IG_ISOTYPE))
  )

# Post filtering missing analysis
missing_analysis_post <- naclus(data_clin_final %>% select(-MPC))

par(mfrow = c(1, 1))
plot(missing_analysis_post)
naplot(missing_analysis_post, which = "all")


#Genomic Dataset Alignement
data_gen_filtered <- data_gen %>%
  semi_join(data_clin_final, by = c("id" = "MPC"))

#constant CNAs are removed
constant_mutations <- data_gen_filtered %>%
  select(-id, -Method) %>%
  select(where(~ n_distinct(., na.rm = TRUE) == 1)) %>% 
  colnames()
print("Zero-variance mutations (constant features):")
print(constant_mutations)
data_gen_filtered <- data_gen_filtered %>%
  select(-all_of(constant_mutations))

# Rare Genomic Mutations (<10%) are removed
print("--- GENOMIC MUTATIONS FREQUENCIES ---")
lapply(data_gen_filtered %>% select(-id, -Method), table)
sapply(data_gen_filtered %>% select(-id, -Method), function(x) {
  mean(x == "1", na.rm = TRUE)
})

rare_mutations <- data_gen_filtered %>%
  select(-id, -Method) %>%
  select(where(~ {
    freq_table <- prop.table(table(.x))
    min(freq_table) < 0.1
  })) %>%
  colnames() %>%
  # Remove del_17p from the exclusion list given its critical high-risk prognostic value.
  setdiff("del_17p")
print("Excluded rare mutations (< 10% prevalence):")
print(rare_mutations)
data_gen_filtered <- data_gen_filtered %>%
  select(-all_of(rare_mutations))
ncol(data_gen_filtered)

#CNAs multicollinearity matrix
cor_matrix <- data_gen_filtered %>%
  select(-id, -Method) %>%
  mutate(across(everything(), as.numeric)) %>%
  cor(use = "pairwise.complete.obs", method = "pearson")

corrplot(cor_matrix, 
         method = "color", 
         type = "upper", 
         order = "hclust", 
         tl.col = "black", 
         tl.cex = 0.6,      
         na.label = "square", 
         na.label.col = "grey")

# Merge the clean clinical and filtered genomic dataset
merged_dataset <- data_clin_final %>%
  left_join(data_gen_filtered %>% select(-Method), by = c("MPC" = "id"))

print(paste("Final dataset dimensions:", nrow(merged_dataset), "patients and", ncol(merged_dataset), "variables"))

# Extract imaging and liquid biopsy variables
data_imaging_to_merge <- data_clin %>%
  select(MPC, any_of(vars_imaging))

# Merge imaging features into the master dataset preserving their missingness (NAs)
merged_dataset_imag <- merged_dataset %>%
  left_join(data_imaging_to_merge, by = "MPC")

print(paste("Final dataset dimensions with imaging:", nrow(merged_dataset_imag), "patients and", ncol(merged_dataset_imag), "variables"))

# Data Export
saveRDS(merged_dataset, 
        file = "Clinical_Genomic_Dataset_Pre_Imputation.rds")
saveRDS(merged_dataset_imag, 
        file = "Imaging_Dataset_Pre_Imputation.rds")


