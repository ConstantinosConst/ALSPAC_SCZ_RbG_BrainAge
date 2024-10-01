###### Project:Assessing the association between global structural brain age and polygenic risk for schizophrenia in early adulthood: A recall-by-genotype study ######
###### Script 2c: ENIGMA brain age model performance check and output preparation for downstream analyses - step 2c
###### Description: Generate metrics/plots for model performance in the current sample and prepares output for downstream analyses
###### Written by: Constantinos Constantinides ######
###### Adapted from: ENIGMA-SCZ-BrainAge study (https://github.com/ConstantinosConst/ENIGMA-SZ-BrainAge.git)
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
            "caret","gridExtra","Hmisc","pastecs","psych","ggplot2", "dplyr","smplot2", "sensemakr"))


# source functions
cat("Prep: sourcing functions\n")

source("get.means.R")  
source("prepare.files.R")
source("model.fits.brainAge.R")
source("model.fits.brainAge.SCZ_PRS.R")

# derive mean volumnes / thickness and surfarea

cat("Reading in data\n")
cat("deriving mean values for thickness, surface and volume\n") 
  
Thick=get.means("CorticalMeasuresENIGMA_ThickAvg.csv") 
Thick$ICV=NULL  

Surf=get.means("CorticalMeasuresENIGMA_SurfAvg.csv") 
Surf$ICV=NULL 

Vol=get.means("SubcorticalMeasuresENIGMA_VolAvg.csv") 

# merged all together
TS=merge(Thick,Surf,by="row.names")

TSV=merge(TS,Vol,by.x="Row.names",by.y="row.names")

# read in covariates
Covs <- read.csv("Covariates.csv"); #Read in the covariates file

# get for missings in SCZ_PRS, Age, Sex
if (length(table(is.na(Covs[,c("SCZ_PRS","Sex","Age")])))>1){
  stop("missing data in SCZ_PRS, Age or Sex")
}

# Check that all of the required columns are present
mcols=c("SubjID","SCZ_PRS","Sex","Age")
colind=match(mcols,names(Covs))
if(length(which(is.na(colind))) > 0){
  stop('At least one of the required columns in your Covariates.csv file is missing. Make sure that the column names are spelled exactly as listed\n
       It is possible that the problem column(s) is: ', mcols[which(is.na(colind))])
}

# Check for duplicated SubjIDs that may cause issues with merging data sets.
if(anyDuplicated(Covs[,c("SubjID")]) != 0) { stop('You have duplicate SubjIDs in your Covariates.csv file.\nMake sure there are no repeat SubjIDs.') }


#combine the files into one dataframe
data = merge(Covs, TSV, by.x="SubjID", by.y="Row.names")

cat("creating csv files for brainAge estimation\n")
df=prepare.files(data,names(Covs))
males=df$males
females=df$females
rm(df)


cat("Step 3: performance metrics for the current sample")

#create log file
sink("ENIGMA_Performance_Metrics.log", type="output")

model.fits.brainAge(males,"males_raw_out.csv")
model.fits.brainAge(females,"females_raw_out.csv")

model.fits.brainAge.SCZ_PRS(males,"males_raw_out.csv")
model.fits.brainAge.SCZ_PRS(females,"females_raw_out.csv")

sink()
## Check that 4 .pdf files and 1 .log file provide performance are stored in your working directory. 

cat("Read-in brain-predicted age for males and females and prepare dataset for downstream analyses")

# males 
output <- read.csv("males_raw_out.csv", header=T, sep="\t")
males$predAge_ENG <- output$age_prediction
rm(output)

# females
output <- read.csv("females_raw_out.csv", header=T, sep="\t")
females$predAge_ENG <- output$age_prediction
rm(output)

BA=rbind(males,females)

# check that no sample duplication
cat("duplicated IDs: ",which(duplicated(BA$SubjID)),"\n")

# merge covariates with brain-predicted age
data=merge(data[,names(Covs),],BA[,c("SubjID","predAge_ENG")], by="SubjID")

# load and merge participants with missing/invalid values on IDPs for descriptive purposes
load("ALSPAC_SCZ_RbG_IDP_Missing.RData")
ALSPAC_SCZ_RbG_IDP_Missing$predAge_ENG <- NA
data <- rbind(data, ALSPAC_SCZ_RbG_IDP_Missing)
rm(ALSPAC_SCZ_RbG_IDP_Missing)

# Save dataset for downstream statistical analyses of PE and brain-PAD (scripts 4a/b)
save(data,file="ALSPAC_SCZ_RbG_BA_ENG.RData")

cat('finished with step 2c!')

### done ### 
### Proceed to step 3a ### 
