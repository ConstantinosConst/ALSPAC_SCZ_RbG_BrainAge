###### Project: Assessing the association between global structural brain age and polygenic risk for schizophrenia in early adulthood: A recall-by-genotype study ######
###### Script 1: Extraction of the ALSPAC-SCZ-RbG imaging sub-sample from the ALSPAC dataset ######
###### Description: Selects participants of the original ALSPAC-SCZ-RbG imaging sub-study along with variables of interest from the wider ALSPAC cohort ######
###### Written by: Constantinos Constantinides ######
###### Last edited: 28/09/2024 ######

# Set working directory to data file location
# setwd("...")

# Load libraries
library(haven)
library(labelled)
library(stringr)
library(dplyr)
library(tidyr)
options(max.print = 1000000)

# Load ALSPAC dataset for approved project (B3067)
df <- read_sav("./CIDB3067_26Aug2023.sav")
## Obs N = 15645, Variables N = 3895

# generate a searchable data dictionary from the full dataset
# ALSPAC_CIDB3067_Dictionary <- labelled::generate_dictionary(df)
## save(ALSPAC_CIDB3067_Dictionary,file="ALSPAC_CIDB3067_26Aug_Dictionary.RData")

# merge the cidB3067 and qlet to derive an ID column for each G1 offspring
df$SubjID <-str_c(df$cidB3067, '', df$qlet)
## cidB3067 provides a unique number per each mother for the approved ALSPAC project (B3067)
## qlet indicates the order of birth (i.e, A, B) of each offspring

# load genetic PCs for approved project
PCs <- read.csv ("B3067_PCs.csv", header=TRUE)
## keep only top 5 PCs
PCs <- PCs[, c("SubjID", "PC1", "PC2", "PC3", "PC4", "PC5")]
# Run a left join to merge with PCs with the core ALSPAC dataset
df  <- merge(x=df, y= PCs, by="SubjID", all.x=TRUE)
rm(PCs)

# Get the total number of G1 offspring participated in the ALSPAC-PE imaging study (MRI.study1)
MRI_Study2_size <- factor(df$MRI.study2)
summary(MRI_Study2_size)
## N=197

# Subset offspring who participated in the ALSPAC-PE imaging sub-study study 
ALSPAC_SCZ_RbG_Imaging <- df[which(df$MRI.study2 == 1),]
rm(df)

# Select variables of interest for the ALSPAC-PE-BraiAge study
## NOTE (SOS): Column names for imaging-derived phenotype (IDPs) may differ in future data releases
## See accompanied excel spreadsheet for column names as currently appear in the ALSPAC dictionary/search tool (https://variables.alspac.bris.ac.uk)

ALSPAC_SCZ_RbG_Imaging <- ALSPAC_SCZ_RbG_Imaging [, c("SubjID", "casecontrol.2", "ageyears.2", "kz021", 
                          "hand", "kz030", "f8ws112", "tc4025a", "FJAL4000", "FJPL163", "FJPL172","CCU3331", "YPB5180", "FKDQ1030", 
                          "FKPL2250", "FKMS1040", "YPF7520", "PC1", "PC2", "PC3", "PC4", "PC5", "MRI_YP_L1157", "MRI_YP_L1158", 
                          "l_bankssts_thickavg.2", "l_caudalanteriorcingulate_thickavg.2", "l_caudalmiddlefrontal_thickavg.2",    
                          "l_cuneus_thickavg.2", "l_entorhinal_thickavg.2", "l_fusiform_thickavg.2", "l_inferiorparietal_thickavg.2",       
                          "l_inferiortemporal_thickavg.2", "l_isthmuscingulate_thickavg.2", "l_lateraloccipital_thickavg.2", "l_lateralorbitofrontal_thickavg.2",   
                          "l_lingual_thickavg.2", "l_medialorbitofrontal_thickavg.2", "l_middletemporal_thickavg.2", "l_parahippocampal_thickavg.2",        
                          "l_paracentral_thickavg.2", "l_parsopercularis_thickavg.2", "l_parsorbitalis_thickavg.2", "l_parstriangularis_thickavg.2",       
                          "l_pericalcarine_thickavg.2", "l_postcentral_thickavg.2", "l_posteriorcingulate_thickavg.2", "l_precentral_thickavg.2",              
                          "l_precuneus_thickavg.2", "l_rostralanteriorcingulate_thick.2", "l_rostralmiddlefrontal_thickavg.2", "l_superiorfrontal_thickavg.2",     
                          "l_superiorparietal_thickavg.2", "l_superiortemporal_thickavg.2", "l_supramarginal_thickavg.2", "l_frontalpole_thickavg.2",            
                          "l_temporalpole_thickavg.2", "l_transversetemporal_thickavg.2", "l_insula_thickavg.2", "r_bankssts_thickavg.2",               
                          "r_caudalanteriorcingulate_thicka.2", "r_caudalmiddlefrontal_thickavg.2", "r_cuneus_thickavg.2", "r_entorhinal_thickavg.2",             
                          "r_fusiform_thickavg.2", "r_inferiorparietal_thickavg.2", "r_inferiortemporal_thickavg.2", "r_isthmuscingulate_thickavg.2",     
                          "r_lateraloccipital_thickavg.2", "r_lateralorbitofrontal_thickavg.2", "r_lingual_thickavg.2", "r_medialorbitofrontal_thickavg.2",  
                          "r_middletemporal_thickavg.2", "r_parahippocampal_thickavg.2", "r_paracentral_thickavg.2", "r_parsopercularis_thickavg.2",       
                          "r_parsorbitalis_thickavg.2", "r_parstriangularis_thickavg.2", "r_pericalcarine_thickavg.2", "r_postcentral_thickavg.2",            
                          "r_posteriorcingulate_thickavg.2", "r_precentral_thickavg.2", "r_precuneus_thickavg.2", "r_rostralanteriorcingulate_thick.2", 
                          "r_rostralmiddlefrontal_thickavg.2", "r_superiorfrontal_thickavg.2", "r_superiorparietal_thickavg.2", "r_superiortemporal_thickavg.2",       
                          "r_supramarginal_thickavg.2","r_frontalpole_thickavg.2", "r_temporalpole_thickavg.2", "r_transversetemporal_thickavg.2", 
                          "r_insula_thickavg.2", "lthickness.2", "rthickness.2", 
                          "l_bankssts_surfavg.2", "l_caudalanteriorcingulate_surfav.2", "l_caudalmiddlefrontal_surfavg.2",   
                          "l_cuneus_surfavg.2", "l_entorhinal_surfavg.2", "l_fusiform_surfavg.2", "l_inferiorparietal_surfavg.2",      
                          "l_inferiortemporal_surfavg.2", "l_isthmuscingulate_surfavg.2", "l_lateraloccipital_surfavg.2", "l_lateralorbitofrontal_surfavg.2", 
                          "l_lingual_surfavg.2", "l_medialorbitofrontal_surfavg.2", "l_middletemporal_surfavg.2", "l_parahippocampal_surfavg.2",      
                          "l_paracentral_surfavg.2", "l_parsopercularis_surfavg.2", "l_parsorbitalis_surfavg.2", "l_parstriangularis_surfavg.2",      
                          "l_pericalcarine_surfavg.2", "l_postcentral_surfavg.2", "l_posteriorcingulate_surfavg.2", "l_precentral_surfavg.2",             
                          "l_precuneus_surfavg.2", "l_rostralanteriorcingulate_surfa.2", "l_rostralmiddlefrontal_surfavg.2", "l_superiorfrontal_surfavg.2",       
                          "l_superiorparietal_surfavg.2", "l_superiortemporal_surfavg.2", "l_supramarginal_surfavg.2", "l_frontalpole_surfavg.2",           
                          "l_temporalpole_surfavg.2", "l_transversetemporal_surfavg.2", "l_insula_surfavg.2", "r_bankssts_surfavg.2",              
                          "r_caudalanteriorcingulate_surfav.2", "r_caudalmiddlefrontal_surfavg.2", "r_cuneus_surfavg.2", "r_entorhinal_surfavg.2",            
                          "r_fusiform_surfavg.2", "r_inferiorparietal_surfavg.2", "r_inferiortemporal_surfavg.2", "r_isthmuscingulate_surfavg.2",      
                          "r_lateraloccipital_surfavg.2", "r_lateralorbitofrontal_surfavg.2", "r_lingual_surfavg.2", "r_medialorbitofrontal_surfavg.2",  
                          "r_middletemporal_surfavg.2", "r_parahippocampal_surfavg.2", "r_paracentral_surfavg.2", "r_parsopercularis_surfavg.2",       
                          "r_parsorbitalis_surfavg.2", "r_parstriangularis_surfavg.2", "r_pericalcarine_surfavg.2", "r_postcentral_surfavg.2",           
                          "r_posteriorcingulate_surfavg.2", "r_precentral_surfavg.2", "r_precuneus_surfavg.2", "r_rostralanteriorcingulate_surfa.2",
                          "r_rostralmiddlefrontal_surfavg.2", "r_superiorfrontal_surfavg.2", "r_superiorparietal_surfavg.2", "r_superiortemporal_surfavg.2",      
                          "r_supramarginal_surfavg.2", "r_frontalpole_surfavg.2", "r_temporalpole_surfavg.2", "r_transversetemporal_surfavg.2",    
                          "r_insula_surfavg.2", "lsurfarea.2", "rsurfarea.2",
                          "llatvent.2", "rlatvent.2","lthal.2", "rthal.2", 
                          "lcaud.2", "rcaud.2", "lput.2", "rput.2", "lpal.2",    
                          "rpal.2", "lhippo.2", "rhippo.2", "lamyg.2", "ramyg.2", 
                          "laccumb.2", "raccumb.2", "icv.2"
                           )] 

# Imaging-derived phenotype were previously derived from sMRI scans by Sharp et al (2020). Read data note for more details: 
# Sharp TH, et al. Population neuroimaging: generation of a comprehensive data resource within the ALSPAC pregnancy and birth cohort. 
# Wellcome Open Res. 2020 Aug 28;5:203. doi: 10.12688/wellcomeopenres.16060.1. PMID: 33043145; PMCID: PMC7531050.


# Rename non-imaging variables/covariates
ALSPAC_SCZ_RbG_Imaging <- ALSPAC_SCZ_RbG_Imaging %>% 
  rename(SCZ_PRS = casecontrol.2,
         Sex = kz021,
         Age = ageyears.2,
         HAND = hand,
         BW = kz030,
         IQ_8 = f8ws112,
         SDQ_ES_17 = tc4025a,
         AUDIT_18 = FJAL4000,
         PE_18_3Level = FJPL163, 
         PsyDis_18 = FJPL172,
         CASTgroup_20 = CCU3331,
         MFQ_22 = YPB5180,
         BMI_24 = FKMS1040, 
         GAD_24 = FKDQ1030,
         PE_24_4level = FKPL2250,
         UniAtt_26 = YPF7520,
         CorticalQC = MRI_YP_L1157,
         SubCorticalQC = MRI_YP_L1158)
colnames(ALSPAC_SCZ_RbG_Imaging)

### Re-code sex
ALSPAC_SCZ_RbG_Imaging$Sex[ALSPAC_SCZ_RbG_Imaging$Sex == 1] <- 0
ALSPAC_SCZ_RbG_Imaging$Sex[ALSPAC_SCZ_RbG_Imaging$Sex == 2] <- 1

# Remove .1 from column names of imaging-derived phenotype (IDPs)
names(ALSPAC_SCZ_RbG_Imaging) <- sub('.2$', '', names(ALSPAC_SCZ_RbG_Imaging))
colnames(ALSPAC_SCZ_RbG_Imaging)
# Fix MFQ_22 and PC2 names
ALSPAC_SCZ_RbG_Imaging <- ALSPAC_SCZ_RbG_Imaging %>% 
  rename(MFQ_22 = MFQ_,
         PC2 = P)
         

# Inspect dataset
summary(ALSPAC_SCZ_RbG_Imaging)

# Remove one participant with 'NA' on SCZ-PRS status (possibly due to withdrawal of consent)
ALSPAC_SCZ_RbG_Imaging <- ALSPAC_SCZ_RbG_Imaging %>% 
  drop_na(SCZ_PRS)

# save dataset for brain age prediction (steps 2 and 3)
write.csv(ALSPAC_SCZ_RbG_Imaging, file="ALSPAC_SCZ_RbG_Imaging.csv", row.names = FALSE)

### End ### 
### Move to script(s) 2a and/or 3a ###