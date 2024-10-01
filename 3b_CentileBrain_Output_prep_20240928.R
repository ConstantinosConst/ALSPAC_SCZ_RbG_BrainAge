###### Project: Assessing the association between global structural brain age and polygenic risk for schizophrenia in early adulthood: A recall-by-genotype study ###### 
###### Script 3b: CentileBrain brain age model performance check and dataset preparation for downstream analyses - step 3b ######
###### Description: Generates performance metrics/plots for the current sample and prepare output for downstream analyses ######
###### Written by: Constantinos Constantinides ######
###### Last edited: 28/09/2024 ######

# Set working directory to file location
# setwd ("...")

# install/load libraries
cat("Prep: installing/loading libraries\n")

load.lib <- function(x){
  for( i in x ){
    if( ! require( i , character.only = TRUE ) ){
      install.packages( i , dependencies = TRUE )
      #require( i , character.only = TRUE )
      library(i)
    }
  }
}

#  Then try/install packages...
load.lib( c("ppcor" , "lsmeans" , "multcomp","data.table","plyr","ModelMetrics",
            "caret","gridExtra","Hmisc","pastecs","psych","ggplot2", "dplyr","smplot2", "sensemakr", "rJava", "xlsx") )


# Source function
source("CB.model.fits.brainAge.R")
source("CB.model.fits.brainAge.SCZ_PRS.R")

# Load ALSPAC-SCZ_RbG imaging dataset (as prepared in script 1)
data <- read.csv("ALSPAC_SCZ_RbG_Imaging.csv", header=TRUE)

# Isolate covariates (non-imaging variables) to merge later with model output
## Non-imaging variables
Covariates <- data[,c("SubjID", "SCZ_PRS", "Age",                               
                      "Sex", "HAND", "BW",                                 
                      "IQ_8", "SDQ_ES_17", "AUDIT_18",                          
                      "PE_18_3Level", "PsyDis_18", "CASTgroup_20",                       
                      "MFQ_22", "GAD_24", "PE_24_4level", "BMI_24", 
                      "UniAtt_26", "PC1", "PC2", "PC3", "PC4", "PC5", 
                      "CorticalQC", "SubCorticalQC" 
)]
### Check for duplicated columns
duplicated(Covariates) 
rm(data)

# Load CB model input files and merge with covariates (as prepared in script 2)
# Males
males <- read.xlsx("CB_males_raw.xlsx", 1)
names(males)[names(males) == "SubjectID"] <- "SubjID"
males <- merge(males, Covariates, by='SubjID')
# females
females <- read.xlsx("CB_females_raw.xlsx", 1)
names(females)[names(females) == "SubjectID"] <- "SubjID"
females <- merge(females, Covariates, by='SubjID')

rm (Covariates)
cat("Step 3: performance metrics for the current sample")

#Create log file for performance metrics
sink("CentileBrain_Performance_Metrics.log", type="output")

model.fits.brainAge(males,"CB_males_raw_out.csv")
model.fits.brainAge(females,"CB_females_raw_out.csv")

model.fits.brainAge.SCZ_PRS(males,"CB_males_raw_out.csv")
model.fits.brainAge.SCZ_PRS(females,"CB_females_raw_out.csv")

sink()

# check that 4 pdf files and 1 log file with performance metrics/plots have been generated in the working directory. 

# Read-in CB output files to get brain-predicted age
# males 
output <- read.csv("CB_males_raw_out.csv", header=T)
males$predAge_CB <- output$x
rm(output)

# females
output <- read.csv("CB_females_raw_out.csv", header=T)
females$predAge_CB <- output$x
rm(output)

BA=rbind(males,females)

# check that no sample duplication
cat("duplicated IDs: ",which(duplicated(BA$SubjID)),"\n")

# Isolate and export covariates and brain-predicted age for downstream statistical analyses
data <- BA %>%  dplyr::select(c("SubjID", "SCZ_PRS", "Age",                               
                        "Sex", "HAND", "BW",                                 
                        "IQ_8", "SDQ_ES_17", "AUDIT_18",                          
                        "PE_18_3Level", "PsyDis_18", "CASTgroup_20",                       
                        "MFQ_22", "GAD_24", "PE_24_4level", "BMI_24", 
                        "UniAtt_26", "PC1", "PC2", "PC3", "PC4", "PC5", 
                        "CorticalQC", "SubCorticalQC", "predAge_CB"))

# Read-in and merge with excluded participants due to missing IDPs
load("ALSPAC_SCZ_RbG_IDP_Missing.RData")
ALSPAC_SCZ_RbG_IDP_Missing$predAge_CB <- NA
data <- rbind(data, ALSPAC_SCZ_RbG_IDP_Missing)

# Export dataset for downstream statistical analyses
save(data,file="ALSPAC_SCZ_RbG_BA_CB.RData")

cat('finished with step 3b')

# Done
# Proceed to step 4a