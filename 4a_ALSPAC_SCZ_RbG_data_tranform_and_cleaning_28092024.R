# Project: Assessing the association between global structural brain age and polygenic risk for schizophrenia in early adulthood: A recall-by-genotype study ###### 
# Script: 4a. Data transformation and cleaning for data visualization and statistical analyses - step 4a
# Description: Re-codes/derives (new) variables and remove observation with failed image QC for downstream analyses
# Script author: Constantinos Constaninides
# Written/last update on: 28/09/2024
# Related paper: https://doi.org/10.1016/j.cortex.2023.11.015

# Load relevant packages
# library (ggplot2)
library(dplyr)
# library(pastecs)
# library(psych)
# library(sensemakr)
# library(olsrr)
# library(smplot2)
# library(effectsize)
# library(tidyverse)

options(max.print=999999)

# setwd to working directory

# Load ALSPAC-PE dataset with ENIGMA-predicted brain age
load('ALSPAC_SCZ_RbG_BA_ENG.RData')
ALSPAC_SCZ_RbG_BA_ENG <- data
rm(data)
# keep SubjID and predicted age only
ALSPAC_SCZ_RbG_BA_ENG <- ALSPAC_SCZ_RbG_BA_ENG[,c("SubjID", "predAge_ENG")]

# Load ALSPAC-PE dataset with CentileBrain-predicted brain age
load("ALSPAC_SCZ_RbG_BA_CB.RData")
ALSPAC_SCZ_RbG_BA_CB <- data
rm(data)

# Merge CentileBrain- and ENIGMA-predicted brain age (plus covariates)
data <- merge(ALSPAC_SCZ_RbG_BA_CB, ALSPAC_SCZ_RbG_BA_ENG, by = "SubjID")

#Inspect dataset
## View(data)
str(data)
head(data)
summary(data)
colnames(data)
## [1] "SubjID"        "SCZ_PRS"       "Age"           "Sex"           "HAND"          "BW"           
## [7] "IQ_8"          "SDQ_ES_17"     "AUDIT_18"      "PE_18_3Level"  "PsyDis_18"     "CASTgroup_20" 
## [13] "MFQ_22"        "GAD_24"        "PE_24_4level"  "BMI_24"        "UniAtt_26"     "PC1"          
## [19] "PC2"           "PC3"           "PC4"           "PC5"           "CorticalQC"    "SubCorticalQC"
## [25] "predAge_CB"    "predAge_ENG"  


# Step 4.1. Fix and derive variables

# Fix coding and/or set as categorical (see accompanying excel file for details on recoding)
## Sex at birth 
data$Sex <- as.factor(data$Sex)
## PEs at age 18
data$SCZ_PRS <- factor(data$SCZ_PRS)
## Handedness
data$HAND[data$HAND == -1] <- NA
data$HAND <- as.factor(data$HAND)
## GAD at age 24
data$GAD_24 <- as.factor(data$GAD_24)
## PE at age 18 (3 level, by adding psychotic disorder as group)
data$PE_18_3Level[data$PE_18_3Level == 1] <- 0
data$PE_18_3Level[data$PE_18_3Level == 2] <- 1
data$PE_18_3Level[data$PE_18_3Level == 3] <- 1
data$PE_18_3Level[data$PsyDis_18 == 1] <- 2
data$PE_18_3Level <- as.factor(data$PE_18_3Level)
data$PsyDis_18 <- NULL
## CAST at age 20 (group; 3-level to 2-level)
data$CASTgroup_20[data$CASTgroup_20 == 1] <- 0
data$CASTgroup_20[data$CASTgroup_20 == 2] <- 0
data$CASTgroup_20[data$CASTgroup_20 == 3] <- 1
data$CASTgroup_20 <- as.factor(data$CASTgroup_20)
## Study for university degree by age 24
data$UniAtt_26[data$UniAtt_26 == 2] <- 1
data$UniAtt_26  <- as.factor(data$UniAtt_26)
## PE at age 24 (4_level to 3-level)
data$PE_24_4level[data$PE_24_4level == 2] <- 1
data$PE_24_4level[data$PE_24_4level == 3] <- 2
names(data)[names(data) == 'PE_24_4level'] <- 'PE_24_3level'
data$PE_24_3level <- as.factor(data$PE_24_3level)
## Image QC variables
data$CorticalQC <- as.factor(data$CorticalQC)
data$SubCorticalQC <- as.factor(data$SubCorticalQC)
## summary(data)

## Derive Brain-PAD (outcome), plus absolute errors

### CentileBrain-derived brain-PAD
data$devAge_CB <- (data$predAge_CB - data$Age)
### Absolute error of CentileBrain brain-age model, Version 2
data$AE_CB <- abs(data$predAge_CB - data$Age)

### ENIGMA_derived brain-PAD
data$devAge_ENG <- (data$predAge_ENG - data$Age)
### Absolute error of ENIGMA brain-age model (for model performance)
data$AE_ENG <- abs(data$predAge_ENG - data$Age)


## Derive new variables (predictors/covariates)

# Derive a cumulative measure of psychotic experiences recorded by age 24 (based on assessments at age 18 and age 24)
data <- data %>% 
  mutate(PE_24_18= case_when(is.na(PE_18_3Level) & is.na(PE_24_3level) ~ NA, 
                             is.na(PE_24_3level) & PE_18_3Level == 0 ~ NA,
                             PE_24_3level == 0 & PE_18_3Level == 0 ~ 0,
                             is.na(PE_18_3Level) & PE_24_3level == 0 ~ 0, 
                             is.na(PE_24_3level) & PE_18_3Level == 1 ~ 1, 
                             is.na(PE_18_3Level) & PE_24_3level == 1 ~ 1,
                             PE_24_3level == 1 & PE_18_3Level == 0 ~ 1,
                             PE_24_3level == 0 & PE_18_3Level == 1 ~ 1,
                             PE_24_3level == 1 & PE_18_3Level == 1 ~ 1,
                             is.na(PE_24_3level) & PE_18_3Level == 2 ~ 2, 
                             is.na(PE_18_3Level) & PE_24_3level == 2 ~ 2,
                             PE_24_3level == 2 & PE_18_3Level == 0 ~ 2,
                             PE_24_3level == 0 & PE_18_3Level == 2 ~ 2,
                             PE_24_3level == 1 & PE_18_3Level == 2 ~ 2,
                             PE_24_3level == 2 & PE_18_3Level == 1 ~ 2,
                             PE_24_3level == 2 & PE_18_3Level == 2 ~ 2))

data$PE_24_18 <- as.factor(data$PE_24_18)
summary(data)


# In addition, convert PE variable from categorical to a binary variable (i.e., with or without PE by age 24) for exploratory analyses
data$PE_24_18_bin <- data$PE_24_18
data$PE_24_18_bin <- as.numeric(data$PE_24_18_bin)
summary(data$PE_24_18_bin)
data$PE_24_18_bin[data$PE_24_18_bin == 1] <- 0
data$PE_24_18_bin[data$PE_24_18_bin == 2] <- 1
data$PE_24_18_bin[data$PE_24_18_bin == 3] <- 1
data$PE_24_18_bin <- as.factor(data$PE_24_18_bin)
summary(data$PE_24_18_bin)


# Step 4.2. Inspection and exclusion of scans that failed image reconstruction/QC

# Select those who failed image reconstruction (= NA) /QC for cortical/sub-cortical QC (=2)
data_exc <- data %>% filter(CorticalQC == 2 | is.na(CorticalQC) | SubCorticalQC == 2 | is.na(SubCorticalQC))
# Get descriptives
summary(data_exc)

## Select those with 'pass' or 'moderate' cortical and subcortical QC for downstream analyses 
data <- data %>% filter((CorticalQC == "0" | CorticalQC == "1") & 
                          (SubCorticalQC == "0" | SubCorticalQC == "1"))
## n=189
# save final dataset
save(data, file='ALSPAC_SCZ_RbG_BA_FINAL.RData')

cat('finished with step 4a')

# done
# Move to script 4b