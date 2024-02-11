##################### Statistical analyses of SCZ-PRS and Brain Age in ALSPAC-SCZ-RBG (Constantinides et al, 2023) #####################
# Author: Constantinos Constaninides
# Written on: 21/04/2023
# Revised on: 15/09/2023
# Last edited (minor corrections) on: 07/01/2024
# Related publication (DOI): https://doi.org/10.1016/j.cortex.2023.11.015

# Load relevant packages
library (ggplot2)
library(dplyr)
library(pastecs)
library(psych)
library(haven)
library(labelled)
library(sensemakr)
library(olsrr)
library(smplot2)

options(max.print=999999)

# Load derived dataset (incl. brain age estimates)
#setwd to working directory
load("ALSPAC_SCZ_RbG_BA_Dataset.RData")
data <- ALSPAC_SCZ_RbG_BA_Dataset
rm(ALSPAC_SCZ_RbG_BA_Dataset)
#Inspect dataset 
str(data)
head(data)
summary(data)
colnames(data)
#remove SPSS label
# var_label(data) <- NULL

# Recode existing variables as needed (see accompanying excel file for details)
data$Sex[data$Sex == 1] <- 0
data$Sex[data$Sex == 2] <- 1
data$HAND[data$HAND == -1] <- NA
data$YPF7520[data$YPF7520 == 2] <- 1
data$CCU3331[data$CCU3331 == 1] <- 0
data$CCU3331[data$CCU3331 == 2] <- 0
data$CCU3331[data$CCU3331 == 3] <- 1
data$FJPL163[data$FJPL163 == 1] <- 0
data$FJPL163[data$FJPL163 == 2] <- 1
data$FJPL163[data$FJPL163 == 3] <- 1
data$FKPL2250[data$FKPL2250== 2] <- 1
data$FKPL2250[data$FKPL2250== 3] <- 2
summary(data)

# Fix variables
### rename and/or convert to categorical 
data$Sex <- as.factor(data$Sex)
data$SCZ_PRS <- as.factor(data$SCZ_PRS)
data$HAND <- as.factor(data$HAND)
data$UniAtt_26 <- as.factor(data$YPF7520)
data$YPF7520 <- NULL
data$CASTgroup_22 <- as.factor(data$CCU3331)
data$CCU3331 <- NULL
data$PE_24 <- as.factor(data$FKPL2250)
data$FKPL2250 <- NULL
data$GAD_24 <- as.factor(data$FKDQ1030)
data$FKDQ1030 <- NULL
data$CorticalQC <- as.factor(data$MRI_YP_L1157)
data$MRI_YP_L1157 <- NULL
data$SubCorticalQC<- as.factor(data$MRI_YP_L1158)
data$MRI_YP_L1158 <- NULL
# merge FJPL163 and FJPL172 into one variable and convert to categorical
data$PE_18 <- data$FJPL163
data$PE_18[data$FJPL172 == 1] <- 2
data$PE_18 <- as.factor(data$PE_18)
## Convert to numerical
data$PC1 <- as.numeric(data$PC1)
data$PC2 <- as.numeric(data$PC2)
data$PC3 <- as.numeric(data$PC3)
data$PC4 <- as.numeric(data$PC4)
data$PC5 <- as.numeric(data$PC5)
#Rename
data$BMI_24 <- data$FKMS1040
data$FKMS1040 <- NULL
data$AUDIT_18 <- data$FJAL4000
data$FJAL4000 <- NULL
data$SDQ_ES_17 <- data$tc4025a
data$tc4025a <- NULL
data$MFQ_22 <- data$YPB5180
data$YPB5180 <- NULL
data$BW <- data$kz030
data$kz030 <-NULL
data$IQ_8 <- data$f8ws112
data$f8ws112 <- NULL

summary(data)
colnames(data)

# Isolate variables needed for the current analyses
data <- data %>% select(SubjID, Sex, Age, SCZ_PRS, predAge_ENG, predAge_CB, PC1, PC2, PC3, PC4, PC5, HAND, UniAtt_26, CASTgroup_22, PE_24, PE_18, GAD_24, CorticalQC, SubCorticalQC, BMI_24, AUDIT_18, SDQ_ES_17, MFQ_22, BW, IQ_8)
summary(data)

# Derive new variables
# ENIGMA_derived brain-PAD
data$devAge_ENG <- (data$predAge_ENG - data$Age)
# absolute error of ENIGMA brain-age model
data$AE_ENG <- abs(data$predAge_ENG - data$Age)
# CentileBrain-derived brain-PAD
data$devAge_CB <- (data$predAge_CB - data$Age)
# Absolute error of CentileBrain brain-age model
data$AE_CB <- abs(data$predAge_CB - data$Age)

# Derive a cumulative measure of psychotic experiences recorded by age 24 (based on assessments at age 18 and age 24)
data <- data %>% 
  mutate(PE_24_18= case_when(is.na(PE_18) & is.na(PE_24) ~ NA, 
                             is.na(PE_24) & PE_18 == 0 ~ NA,
                             PE_24 == 0 & PE_18 == 0 ~ 0,
                             is.na(PE_18) & PE_24 == 0 ~ 0, 
                             is.na(PE_24) & PE_18 == 1 ~ 1, 
                             is.na(PE_18) & PE_24 == 1 ~ 1,
                             PE_24 == 1 & PE_18 == 0 ~ 1,
                             PE_24 == 0 & PE_18 == 1 ~ 1,
                             PE_24 == 1 & PE_18 == 1 ~ 1,
                             is.na(PE_24) & PE_18 == 2 ~ 2, 
                             is.na(PE_18) & PE_24 == 2 ~ 2,
                             PE_24 == 2 & PE_18 == 0 ~ 2,
                             PE_24 == 0 & PE_18 == 2 ~ 2,
                             PE_24 == 1 & PE_18 == 2 ~ 2,
                             PE_24 == 2 & PE_18 == 1 ~ 2,
                             PE_24 == 2 & PE_18 == 2 ~ 2))

data$PE_24_18 <- as.factor(data$PE_24_18)
summary(data)

# In addition, convert PE variable from categorical to a binary variable (i.e., with or without PE by age 24) for explatory analyses
data$PE_24_18_bin <- data$PE_24_18
data$PE_24_18_bin <- as.numeric(data$PE_24_18_bin)
summary(data$PE_24_18_bin)
data$PE_24_18_bin[data$PE_24_18_bin == 1] <- 0
data$PE_24_18_bin[data$PE_24_18_bin == 2] <- 1
data$PE_24_18_bin[data$PE_24_18_bin == 3] <- 1
data$PE_24_18_bin <- as.factor(data$PE_24_18_bin)
summary(data$PE_24_18_bin)

# Excluded those who failed image QC
data <- data %>% filter(CorticalQC == 0 | CorticalQC == 1)
data <- data %>% filter(SubCorticalQC == 0 | SubCorticalQC == 1)

# Get summary statistics per PRS group for Table 1 and brain-PAD
data_low <- data[which(data$SCZ_PRS == 0),]
data_high <- data[which(data$SCZ_PRS == 1),]
summary(data_low)
summary(data_high)
# Get more stats for continues variables
describe(data_low)
describe(data_high)
stat.desc(data_low)
stat.desc(data_high)

# Run between-group comparisons for Table 1
# Continues variables 
wilcox.test(data$Age ~ data$SCZ_PRS, alternative = "two.sided")
wilcox.test(data$BW ~ data$SCZ_PRS, alternative = "two.sided")
wilcox.test(data$BMI_24 ~ data$SCZ_PRS, alternative = "two.sided")
wilcox.test(data$IQ_8 ~ data$SCZ_PRS, alternative = "two.sided")
wilcox.test(data$SDQ_ES_17 ~ data$SCZ_PRS, alternative = "two.sided")
wilcox.test(data$MFQ_22 ~ data$SCZ_PRS, alternative = "two.sided")
wilcox.test(data$AUDIT_18 ~ data$SCZ_PRS, alternative = "two.sided")
# Categorical
chisq.test(data$Sex, data$SCZ_PRS, correct=FALSE)
chisq.test(data$HAND, data$SCZ_PRS, correct=FALSE)
chisq.test(data$UniAtt_26, data$SCZ_PRS, correct=FALSE)
fisher.test(data$GAD_24, data$SCZ_PRS)
fisher.test(data$PE_24_18, data$SCZ_PRS)
fisher.test(data$CASTgroup_22, data$SCZ_PRS)

# Calculate performance metrics for STAble B1
## With respect to sex
data_M <- data[which(data$Sex == 0),]
data_F <- data[which(data$Sex == 1),]
### MAE/RSD - Males
stat.desc(data_M$AE_ENG)
stat.desc(data_M$AE_CB)
### MAE/RSD - Females
stat.desc(data_F$AE_ENG)
stat.desc(data_F$AE_CB)
## Pearson's R and R2 - Males
cor(data_M$Age,data_M$predAge_ENG)
cor(data_M$Age,data_M$predAge_CB)
caret::R2(data_M$Age,data_M$predAge_ENG)
caret::R2(data_M$Age,data_M$predAge_CB)
## Pearson's R and R2 - Females
cor(data_F$Age,data_F$predAge_ENG)
cor(data_F$Age,data_F$predAge_CB)
caret::R2(data_F$Age,data_F$predAge_ENG)
caret::R2(data_F$Age,data_F$predAge_CB)

## With respect to both sex and SCZ_PRS
data_M_low <- data_M[which(data_M$SCZ_PRS == 0),]
data_M_high <- data_M[which(data_M$SCZ_PRS == 1),]
data_F_low <- data_F[which(data_F$SCZ_PRS == 0),]
data_F_high <- data_F[which(data_F$SCZ_PRS == 1),]
## Low PRS
### MAE/RSD - Males / low PRS
stat.desc(data_M_low$AE_ENG)
stat.desc(data_M_low$AE_CB)
### MAE/RSD - Females / low PRS
stat.desc(data_F_low$AE_ENG)
stat.desc(data_F_low$AE_CB)
## Pearson's R and R2 - Males / low PRS
cor(data_M_low$Age,data_M_low$predAge_ENG)
cor(data_M_low$Age,data_M_low$predAge_CB)
caret::R2(data_M_low$Age,data_M_low$predAge_ENG)
caret::R2(data_M_low$Age,data_M_low$predAge_CB)
## Pearson's R and R2 - Females / low PRS
cor(data_F_low$Age,data_F_low$predAge_ENG)
cor(data_F_low$Age,data_F_low$predAge_CB)
caret::R2(data_F_low$Age,data_F_low$predAge_ENG)
caret::R2(data_F_low$Age,data_F_low$predAge_CB)

# High PRS
### MAE/RSD - Males / high PRS
stat.desc(data_M_high$AE_ENG)
stat.desc(data_M_high$AE_CB)
### MAE/RSD - Females / high PRS
stat.desc(data_F_high$AE_ENG)
stat.desc(data_F_high$AE_CB)
## Pearson's R and R2 - Males / high PRS
cor(data_M_high$Age,data_M_high$predAge_ENG)
cor(data_M_high$Age,data_M_high$predAge_CB)
caret::R2(data_M_high$Age,data_M_high$predAge_ENG)
caret::R2(data_M_high$Age,data_M_high$predAge_CB)
## Pearson's R and R2 - Females / low PRS
cor(data_F_high$Age,data_F_high$predAge_ENG)
cor(data_F_high$Age,data_F_high$predAge_CB)
caret::R2(data_F_high$Age,data_F_high$predAge_ENG)
caret::R2(data_F_high$Age,data_F_high$predAge_CB)

# Plot predicted age for each brain age model and with respect to sex (SFigure B1)
colnames(data)
#data_2 <- data[,c(1,2,4,11,12)]
data_2 <- data %>% select(SubjID, Sex, Age, predAge_ENG, predAge_CB)
data_2$Age_type <- ""
data_2_M <- data_2[which(data_2$Sex == 0),]
data_2_F <- data_2[which(data_2$Sex == 1),]
colnames(data_2_M)
data_M_Age <- data_2_M %>% select(SubjID, Sex, Age, Age_type)
data_M_Age$Age_type[data_M_Age$Age_type==""] <- "Chronological age (Males)"
data_M_predAgeENG <- data_2_M %>% select(SubjID, Sex, predAge_ENG, Age_type)
names(data_M_predAgeENG)[names(data_M_predAgeENG) == "predAge_ENG"] <- "Age"
data_M_predAgeENG$Age_type[data_M_predAgeENG$Age_type==""] <- "ENIGMA model-predicted age (Males)"
data_M_predAgeCB <- data_2_M %>% select(SubjID, Sex, predAge_CB, Age_type)
names(data_M_predAgeCB)[names(data_M_predAgeCB) == "predAge_CB"] <- "Age"
data_M_predAgeCB$Age_type[data_M_predAgeCB$Age_type==""] <- "CentileBrain model-predicted age (Males)"
colnames(data_2_F)
data_F_Age <- data_2_F%>% select(SubjID, Sex, Age, Age_type)
data_F_Age$Age_type[data_F_Age$Age_type==""] <- "Chronological age (Females)"
data_F_predAgeENG <- data_2_F %>% select(SubjID, Sex, predAge_ENG, Age_type)
names(data_F_predAgeENG)[names(data_F_predAgeENG) == "predAge_ENG"] <- "Age"
data_F_predAgeENG$Age_type[data_F_predAgeENG$Age_type==""] <- "ENIGMA model-predicted age (Females)"
data_F_predAgeCB <- data_2_F %>% select(SubjID, Sex, predAge_CB, Age_type)
names(data_F_predAgeCB)[names(data_F_predAgeCB) == "predAge_CB"] <- "Age"
data_F_predAgeCB$Age_type[data_F_predAgeCB$Age_type==""] <- "CentileBrain model-predicted age (Females)"

data_plotB1 <- rbind(data_M_Age, data_M_predAgeENG)
data_plotB1 <- rbind(data_plotB1, data_M_predAgeCB)
data_plotB1 <- rbind(data_plotB1, data_F_Age)
data_plotB1 <- rbind(data_plotB1, data_F_predAgeENG)
data_plotB1 <- rbind(data_plotB1, data_F_predAgeCB)

ggplot(data_plotB1, aes(x = Age, color = Age_type, fill = Age_type)) + 
  geom_density(alpha = 0.2) +
  theme_bw() + theme(axis.text=element_text(size=12),
                     axis.title=element_text(size=12),legend.title=element_blank(),
                     legend.text = element_text(size =10),legend.position="bottom") 

# Summary statics for SFigure B1
# Males
describe(data_2_M)
# Females
describe(data_2_F)

# Assess presence of age-related bias in brain-PAD (SFigure 2 A and B)
# ENIGMA-derived brain-PAD
ggplot(data, aes(x=Age, y=devAge_ENG, color=as.factor(SCZ_PRS))) + 
  geom_hline(yintercept=0 ,linetype=2) +
  geom_point(alpha=0.5, size=2) + 
  scale_y_continuous(limits=(c(-20,25))) +
  geom_smooth(method=lm, colour="black", size=0.5, se=FALSE, fullrange=TRUE) +
  geom_smooth(aes(group=SCZ_PRS), method=lm, size=0.5, se=FALSE, fullrange=TRUE, linetype="dashed") +
  scale_color_manual(labels = c("Low PRS", "High PRS"), values = c("#009E73", "#D55E00"), name = "SCZ-PRS") +
  labs(x = "Chronological age", y = "ENIGMA-derived brain-PAD") +
  theme_bw()
cor.test(data$devAge_ENG, data$Age)

# CentileBrain-derived brain-PAD
ggplot(data, aes(x=Age, y=devAge_CB, color=as.factor(SCZ_PRS))) + 
  geom_point(alpha=0.5, size=2) + 
  geom_hline(yintercept=0 ,linetype=2) +
  scale_y_continuous(limits=(c(-3,3))) +
  geom_smooth(method=lm, colour="black", size=0.5, se=FALSE, fullrange=TRUE) +
  geom_smooth(aes(group=SCZ_PRS), method=lm, size=0.5, se=FALSE, fullrange=TRUE, linetype="dashed") +
  scale_color_manual(labels = c("Low PRS", "High PRS"), values = c("#009E73", "#D55E00"), name = "SCZ-PRS") +
  labs(x = "Chronological age", y = "CentileBrain-derived brain-PAD") +
  theme_bw()
cor.test(data$devAge_CB, data$Age)

# residualize Brain-PAD for age to illustrate age-bias correction (SFig. B2, C and D)

# ENIGMA-derived brain-PAD residualised for age 
devAge_ENG_Age = lm(devAge_ENG ~ Age, data=data)
data$devAge_ENG_resAge = resid(devAge_ENG_Age)
# Plot (SFig. B2, C)
ggplot(data, aes(x=Age, y=devAge_ENG_resAge, color=as.factor(SCZ_PRS))) + 
  geom_hline(yintercept=0 ,linetype=2) +
  geom_point(alpha=0.5, size=2) + 
  scale_y_continuous(limits=(c(-20,25))) +
  geom_smooth(method=lm, colour="black", size=0.5, se=FALSE, fullrange=TRUE) +
  geom_smooth(aes(group=SCZ_PRS), method=lm, size=0.5, se=FALSE, fullrange=TRUE, linetype="dashed") +
  scale_color_manual(labels = c("Low PRS", "High PRS"), values = c("#009E73", "#D55E00"), name = "SCZ-PRS") +
  labs(x = "Chronological age", y = "ENIGMA-derived brain-PAD resid. for age") +
  theme_bw()
cor.test(data$devAge_ENG_resAge, data$Age)


# CentileBrain-derived brain-PAD residualised for age 
devAge_CB_Age = lm(devAge_CB ~ Age, data=data)
data$devAge_CB_resAge=resid(devAge_CB_Age)
# Plot (SFig. B2, D)
ggplot(data, aes(x=Age, y=devAge_CB_resAge, color=as.factor(SCZ_PRS))) + 
  geom_hline(yintercept=0 ,linetype=2) +
  geom_point(alpha=0.5, size=2) + 
  scale_y_continuous(limits=(c(-3,3))) +
  geom_smooth(method=lm, colour="black", size=0.5, se=FALSE, fullrange=TRUE) +
  geom_smooth(aes(group=SCZ_PRS), method=lm, size=0.5, se=FALSE, fullrange=TRUE, linetype="dashed") +
  scale_color_manual(labels = c("Low PRS", "High PRS"), values = c("#009E73", "#D55E00"), name = "SCZ-PRS") +
  labs(x = "Chronological age", y = "CentileBrain-derived brain-PAD resid. for age") +
  theme_bw()
cor.test(data$devAge_CB_resAge, data$Age)

# Assess correlation between ENIGMA-derived brain-PAD and Centile-derived brain-PAD (SFig. B3, left)
ggplot(data, aes(x=devAge_CB,y=devAge_ENG, color=as.factor(SCZ_PRS))) + 
  geom_point(alpha=0.5, size=2) + 
  geom_smooth(method=lm, colour="black", size=0.5, se=FALSE, fullrange=TRUE) +
  scale_color_manual(labels = c("Low PRS", "High PRS"), values = c("#009E73", "#D55E00"), name = "SCZ-PRS") +
  labs(x = "CentileBrain-derived brain-PAD", y = "ENIGMA-derived brain-PAD") +
  theme_bw()
cor.test (data$devAge_ENG, data$devAge_CB)

# repeated with age-residualised brain-PAD (SFig. B3, right)
ggplot(data, aes(x=devAge_CB_resAge, y=devAge_ENG_resAge, color=as.factor(SCZ_PRS))) + 
  geom_point(alpha=0.5, size=2) + 
  geom_smooth(method=lm, colour="black", size=0.5, se=FALSE, fullrange=TRUE) +
  scale_color_manual(labels = c("Low PRS", "High PRS"), values = c("#009E73", "#D55E00"), name = "SCZ-PRS") +
  labs(x = "CentileBrain-derived brain-PAD resid. for age", y = "ENIGMA-derived brain-PAD resid. for age") +
  theme_bw()
cor.test (data$devAge_CB_resAge, data$devAge_ENG_resAge)

# Run primary models for difference in brain-PAD between low and high SCZ-PRS (primary analyses / STable B2)

# ENIGMA-dervived brain-PAD
DevAge_SCZ_PRS_Age_Sex=lm(devAge_ENG ~ as.factor(SCZ_PRS) + as.factor(Sex) + Age, data = data)
summary(DevAge_SCZ_PRS_Age_Sex)
confint(DevAge_SCZ_PRS_Age_Sex)
partial_r2(DevAge_SCZ_PRS_Age_Sex)
# Cohen's d
# d<-t.val*(n1+n2)/(sqrt(n1*n2)*sqrt(df))
d_ENG <-(-0.231)*(93+96)/(sqrt(93*96)*sqrt(185))
d_ENG
# check assumptions of linear regression
ols_plot_resid_qq(DevAge_SCZ_PRS_Age_Sex)
ols_test_normality(DevAge_SCZ_PRS_Age_Sex)
ols_test_correlation(DevAge_SCZ_PRS_Age_Sex)
ols_plot_resid_fit(DevAge_SCZ_PRS_Age_Sex)
ols_plot_resid_hist(DevAge_SCZ_PRS_Age_Sex)

# Raincloud plot for ENIGMA-derived brain_PAD in low and high SCZ-PRS, residualised for age and sex (Figure 1)
devAge_ENG_Age_Sex = lm(devAge_ENG ~ Age + Sex, data=data)
data$devAge_ENG_resAgeSex = resid(devAge_ENG_Age_Sex)

Group <- factor(data$SCZ_PRS, labels = c("Low SCZ-PRS", "High SCZ-PRS"))
ggplot(data = data, mapping = aes(x=Group, y=devAge_ENG_resAgeSex, fill= Group)) +
  sm_raincloud(sep_level = 1, point.params = list(size = 3, shape = 21,
                                                  color = 'transparent', alpha = 0.4)) +
  labs(x = "", y = "Brain-PAD resid. for age and sex") +
  scale_fill_manual(labels = c("Low SCZ-PRS", "High SCZ-PRS"), values = c("#009E73","#D55E00"), name = "SCZ-PRS") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(size=12),
                     axis.title=element_text(size=12), legend.title = element_text(size = 12),
                     legend.text = element_text(size = 12), legend.position="none")

#CentileBrain-brain-PAD
DevAge_SCZ_PRS_Age_Sex_CB=lm(devAge_CB ~ as.factor(SCZ_PRS) + as.factor(Sex) + Age, data = data)
summary(DevAge_SCZ_PRS_Age_Sex_CB)
confint(DevAge_SCZ_PRS_Age_Sex_CB)
partial_r2(DevAge_SCZ_PRS_Age_Sex_CB)
# Cohen's d
# d<-t.val*(n1+n2)/(sqrt(n1*n2)*sqrt(df))
#CB
d_CB <-0.195*(93+96)/(sqrt(93*96)*sqrt(185))
d_CB
# check assumptions of linear regression
ols_plot_resid_qq(DevAge_SCZ_PRS_Age_Sex_CB)
ols_test_normality(DevAge_SCZ_PRS_Age_Sex_CB)
ols_test_correlation(DevAge_SCZ_PRS_Age_Sex_CB)
ols_plot_resid_fit(DevAge_SCZ_PRS_Age_Sex_CB)
ols_plot_resid_hist(DevAge_SCZ_PRS_Age_Sex_CB)

## Raincloud plot for CentileBrain-derived brain_PAD in low and high SCZ-PRS, residualised for age and sex (SFigure B4)
devAge_CB_Age_Sex = lm(devAge_CB ~ Age + Sex, data=data)
data$devAge_CB_resAgeSex = resid(devAge_CB_Age_Sex)

ggplot(data = data, mapping = aes(x=Group, y=devAge_CB_resAgeSex, fill= Group)) +
  sm_raincloud(sep_level = 1, point.params = list(size = 3, shape = 21,
                                                  color = 'transparent', alpha = 0.4)) +
  labs(x = "", y = "Brain-PAD resid. for age and sex") +
  scale_fill_manual(labels = c("Low SCZ-PRS", "High SCZ-PRS"), values = c("#009E73","#D55E00"), name = "SCZ-PRS") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(size=12),
                     axis.title=element_text(size=12), legend.title = element_text(size = 12),
                     legend.text = element_text(size = 12), legend.position="none")

# Sensitivity analyses for genetic PCs and brain-PAD outliers

## 1. Addition of top 5 PCs as covariates

###ENIGMA-derived brain-PAD
DevAge_SCZ_PRS_Age_Sex_PCs=lm(devAge_ENG ~ as.factor(SCZ_PRS) + as.factor(Sex) + Age + PC1 + PC2 + PC3 + PC4 + PC5, data = data)
summary(DevAge_SCZ_PRS_Age_Sex_PCs)
confint(DevAge_SCZ_PRS_Age_Sex_PCs)
partial_r2(DevAge_SCZ_PRS_Age_Sex_PCs)
d_ENG_PCs <-(-0.094)*(93+96)/(sqrt(93*96)*sqrt(171))
d_ENG_PCs
#### check assumptions of linear regression
ols_plot_resid_qq(DevAge_SCZ_PRS_Age_Sex_PCs)
ols_test_normality(DevAge_SCZ_PRS_Age_Sex_PCs)
ols_test_correlation(DevAge_SCZ_PRS_Age_Sex_PCs)
ols_plot_resid_fit(DevAge_SCZ_PRS_Age_Sex_PCs)
ols_plot_resid_hist(DevAge_SCZ_PRS_Age_Sex_PCs)

### CB
DevAge_SCZ_PRS_Age_Sex_PCs_CB=lm(devAge_CB ~ as.factor(SCZ_PRS) + as.factor(Sex) + Age + PC1 + PC2 + PC3 + PC4 + PC5, data = data)
summary(DevAge_SCZ_PRS_Age_Sex_PCs_CB)
confint(DevAge_SCZ_PRS_Age_Sex_PCs_CB)
partial_r2(DevAge_SCZ_PRS_Age_Sex_PCs_CB)
d_CB_PCs <- 0.486*(93+96)/(sqrt(93*96)*sqrt(171))
d_CB_PCs

#### check assumptions of linear regression
ols_plot_resid_qq(DevAge_SCZ_PRS_Age_Sex_PCs_CB)
ols_test_normality(DevAge_SCZ_PRS_Age_Sex_PCs_CB)
ols_test_correlation(DevAge_SCZ_PRS_Age_Sex_PCs_CB)
ols_plot_resid_fit(DevAge_SCZ_PRS_Age_Sex_PCs_CB)
ols_plot_resid_hist(DevAge_SCZ_PRS_Age_Sex_PCs_CB)


# 2. Exclusion of outlier for ENIGMA-derived brain-PAD

#Check for outliers in ENIGMA-Brain-PAD (i.e., above or below 3.00 SD from the mean per SCZ-PRS group)
# plot distribution of brian-PAD across the high SCZ-PRS group
ggplot(data_high, aes(x=devAge_ENG)) + 
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.2, fill="blue") +
  labs(x = "ENIGMA-derived brain-PAD")
# Now, Standardize brain-PAD across the high SCZ-PRS group
data_high$devAge_std_ENG <- data_high$devAge_ENG
data_high <- data_high%>% mutate_at(c('devAge_std_ENG'), ~(scale(.) %>% as.vector))
# Plot standardized brain-PAD across the high PRS group
ggplot(data_high, aes(x=devAge_std_ENG)) + 
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.2, fill="blue") +
  labs(x = "standardised ENIGMA-derived brain-PAD")
# Exclude one participant with brain-PAD < -3.00 SDs
data_high_wout <- data_high %>% filter(devAge_std_ENG > (-3.00))
# repeat for low PRS group
ggplot(data_low, aes(x=devAge_ENG)) + 
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.2, fill="blue") +
  labs(x = "ENIGMA-derived brain-PAD")
# Now, Standardize brain-PAD across the high SCZ-PRS group
data_low$devAge_std_ENG <- data_low$devAge_ENG
data_low <- data_low%>% mutate_at(c('devAge_std_ENG'), ~(scale(.) %>% as.vector))
# Plot standardized brain-PAD across the high PRS group
ggplot(data_low, aes(x=devAge_std_ENG)) + 
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.2, fill="blue") +
  labs(x = "standardised ENIGMA-derived brain-PAD")
# No outlier in the low SCZ-PRS group

# Get mean(SD) ENIGMA-derived brain-PAD after excluding one outlier
#Low PRS
describe(data_low$devAge_ENG)
#High PRS
describe(data_high_wout$devAge_ENG)

# Run primary model after exclusion of one outlier in the high-PRS group (N=188)
data_wout <- rbind(data_low, data_high_wout)
DevAge_SCZ_PRS_Age_Sex_wout=lm(devAge_ENG ~ as.factor(SCZ_PRS) + as.factor(Sex) + Age, data = data_wout)
summary(DevAge_SCZ_PRS_Age_Sex_wout)
confint(DevAge_SCZ_PRS_Age_Sex_wout)
partial_r2(DevAge_SCZ_PRS_Age_Sex_wout)
d_wout <- (-0.018)*(93+95)/(sqrt(93*95)*sqrt(18))
d_wout


# Run secondary model addition adjusting for genetic PCs and after exclusion of outlier in the high-PRS group (N=179)
DevAge_SCZ_PRS_Age_Sex_PCs_wout=lm(devAge_ENG ~ as.factor(SCZ_PRS) + as.factor(Sex) + Age + PC1 + PC2 + PC3 + PC4 + PC5, data = data_wout)
summary(DevAge_SCZ_PRS_Age_Sex_PCs_wout)
confint(DevAge_SCZ_PRS_Age_Sex_PCs_wout)
partial_r2(DevAge_SCZ_PRS_Age_Sex_PCs_wout)
d_wout_PCs <- 0.07*(93+95)/(sqrt(93*95)*sqrt(184))
d_wout_PCs

## Check for outlier in CentileBrain-derived brain-PAD
### plot distribution of brian-PAD across the high SCZ-PRS group
ggplot(data_high, aes(x=devAge_CB)) + 
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.2, fill="blue") +
  labs(x = "CentileBrain-derived brain-PAD")
### Now, Standardize brain-PAD across the high SCZ-PRS group
data_high$devAge_std_CB <- data_high$devAge_CB
data_high <- data_high%>% mutate_at(c('devAge_std_CB'), ~(scale(.) %>% as.vector))
### Plot standardized brain-PAD across the high PRS group
ggplot(data_high, aes(x=devAge_std_CB)) + 
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.2, fill="blue") +
  labs(x = "standardised CentileBrain-derived brain-PAD")
### No outliers in the high SCZ-PRS group
### repeat for low PRS group
ggplot(data_low, aes(x=devAge_CB)) + 
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.2, fill="blue") +
  labs(x = "CentileBrain-derived brain-PAD")
### Now, Standardize brain-PAD across the low SCZ-PRS group
data_low$devAge_std_CB <- data_low$devAge_CB
data_low <- data_low%>% mutate_at(c('devAge_std_CB'), ~(scale(.) %>% as.vector))
### Plot standardized brain-PAD across the high PRS group
ggplot(data_low, aes(x=devAge_std_CB)) + 
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.2, fill="blue") +
  labs(x = "standardised CentileBrain-derived brain-PAD")
### No outlier in the low SCZ-PRS group


# Exploratory analyses for the effects of phenotype of interest (IV) and their interactions with SCZ-PRS (STable B3)

#IV: birth weight
### DV: ENIGMA-derived brain-PAD
DevAge_ENG_BW=lm(devAge_ENG ~  Age + as.factor(Sex) + BW + as.factor(SCZ_PRS) + as.factor(SCZ_PRS)*BW, data = data)
summary(DevAge_ENG_BW)
# check assumptions of linear regression
ols_plot_resid_qq(DevAge_ENG_BW)
ols_test_normality(DevAge_ENG_BW)
ols_test_correlation(DevAge_ENG_BW)
ols_plot_resid_fit(DevAge_ENG_BW)
ols_plot_resid_hist(DevAge_ENG_BW)

### DV: CentileBrain-derived brain-PAD
DevAge_CB_BW=lm(devAge_CB ~ Age + as.factor(Sex) + BW + as.factor(SCZ_PRS) + as.factor(SCZ_PRS)*BW, data = data)
summary(DevAge_CB_BW)
# check assumptions of linear regression
ols_plot_resid_qq(DevAge_CB_BW)
ols_test_normality(DevAge_CB_BW)
ols_test_correlation(DevAge_CB_BW)
ols_plot_resid_fit(DevAge_CB_BW)
ols_plot_resid_hist(DevAge_CB_BW)

## IV: IQ at age 8
### DV: ENIGMA-dervived brain-PAD
DevAge_ENG_IQ=lm(devAge_ENG ~  Age + as.factor(Sex) + IQ_8 + as.factor(SCZ_PRS) + as.factor(SCZ_PRS)*IQ_8, data = data)
summary(DevAge_ENG_IQ)
# check assumptions of linear regression
ols_plot_resid_qq(DevAge_ENG_IQ)
ols_test_normality(DevAge_ENG_IQ)
ols_test_correlation(DevAge_ENG_IQ)
ols_plot_resid_fit(DevAge_ENG_IQ)
ols_plot_resid_hist(DevAge_ENG_IQ)

### DV: CentileBrain-derived brain-PAD
DevAge_CB_IQ=lm(devAge_CB ~  Age + as.factor(Sex) + IQ_8 + as.factor(SCZ_PRS) + as.factor(SCZ_PRS)*IQ_8, data = data)
summary(DevAge_CB_IQ)
# check assumptions of linear regression
ols_plot_resid_qq(DevAge_CB_IQ)
ols_test_normality(DevAge_CB_IQ)
ols_test_correlation(DevAge_CB_IQ)
ols_plot_resid_fit(DevAge_CB_IQ)
ols_plot_resid_hist(DevAge_CB_IQ)

# IV: SDQ emotional symptoms score at age 17
## DV: CentileBrain-derived brain-PAD
DevAge_ENG_SDQ=lm(devAge_ENG ~ Age + as.factor(Sex) + SDQ_ES_17+ as.factor(SCZ_PRS) + as.factor(SCZ_PRS)*SDQ_ES_17, data = data)
summary(DevAge_ENG_SDQ)
confint(DevAge_ENG_SDQ)
# check assumptions of linear regression
ols_plot_resid_qq(DevAge_ENG_SDQ)
ols_test_normality(DevAge_ENG_SDQ)
ols_test_correlation(DevAge_ENG_SDQ)
ols_plot_resid_fit(DevAge_ENG_SDQ)
ols_plot_resid_hist(DevAge_ENG_SDQ)

## DV: CentileBrain-derived brain-PAD
DevAge_CB_SDQ=lm(devAge_CB ~ Age + as.factor(Sex) + SDQ_ES_17 + as.factor(SCZ_PRS) + as.factor(SCZ_PRS)*SDQ_ES_17, data = data)
summary(DevAge_CB_SDQ)
confint(DevAge_CB_SDQ)
### check assumption of linear regression
ols_plot_resid_qq(DevAge_CB_SDQ)
ols_test_normality(DevAge_CB_SDQ)
ols_test_correlation(DevAge_CB_SDQ)
ols_plot_resid_fit(DevAge_CB_SDQ)
ols_plot_resid_hist(DevAge_CB_SDQ)

# Depression at age 22
### DV: ENIGMA-derived brain-PAD
DevAge_ENG_MFQ=lm(devAge_ENG ~ Age + Sex + MFQ_22 + as.factor(SCZ_PRS) + as.factor(SCZ_PRS)*MFQ_22 , data = data)
summary(DevAge_ENG_MFQ)
### check assumptions of linear regression
ols_plot_resid_qq(DevAge_ENG_MFQ)
ols_test_normality(DevAge_ENG_MFQ)
ols_test_correlation(DevAge_ENG_MFQ)
ols_plot_resid_fit(DevAge_ENG_MFQ)
ols_plot_resid_hist(DevAge_ENG_MFQ)

### DV: CentileBrain-dervied brain-PAD
DevAge_CB_MFQ=lm(devAge_CB ~ Age + Sex + MFQ_22 + as.factor(SCZ_PRS) + as.factor(SCZ_PRS)*MFQ_22, data = data)
summary(DevAge_CB_MFQ)
### check assumptions of linear regression
ols_plot_resid_qq(DevAge_CB_MFQ)
ols_test_normality(DevAge_CB_MFQ)
ols_test_correlation(DevAge_CB_MFQ)
ols_plot_resid_fit(DevAge_CB_MFQ)
ols_plot_resid_hist(DevAge_CB_MFQ)
                   

# AUDIT at age 18
## ENIGMA-derived brain-PAD
DevAge_ENG_AUDIT=lm(devAge_ENG ~ Age + as.factor(Sex) + AUDIT_18 + as.factor(SCZ_PRS) + as.factor(SCZ_PRS)*AUDIT_18, data = data)
summary(DevAge_ENG_AUDIT)
### check assumptions of linear regression
ols_plot_resid_qq(DevAge_ENG_AUDIT)
ols_test_normality(DevAge_ENG_AUDIT)
ols_test_correlation(DevAge_ENG_AUDIT)
ols_plot_resid_fit(DevAge_ENG_AUDIT)
ols_plot_resid_hist(DevAge_ENG_AUDIT)

DevAge_CB_AUDIT=lm(devAge_CB ~ Age + as.factor(Sex) + AUDIT_18 + as.factor(SCZ_PRS) + as.factor(SCZ_PRS)*AUDIT_18, data = data)
summary(DevAge_CB_AUDIT)
### check assumptions of linear regression
ols_plot_resid_qq(DevAge_CB_AUDIT)
ols_test_normality(DevAge_CB_AUDIT)
ols_test_correlation(DevAge_CB_AUDIT)
ols_plot_resid_fit(DevAge_CB_AUDIT)
ols_plot_resid_hist(DevAge_CB_AUDIT)

## IV: BMI at age 24
### DV: ENIGMA-dervived brain-PAD
DevAge_ENG_BMI=lm(devAge_ENG ~  Age + as.factor(Sex) + BMI_24 + as.factor(SCZ_PRS) + as.factor(SCZ_PRS)*BMI_24, data = data)
summary(DevAge_ENG_BMI)
### check assumptions of linear regression
ols_plot_resid_qq(DevAge_ENG_BMI)
ols_test_normality(DevAge_ENG_BMI)
ols_test_correlation(DevAge_ENG_BMI)
ols_plot_resid_fit(DevAge_ENG_BMI)
ols_plot_resid_hist(DevAge_ENG_BMI)

### DV: CentileBrain-dervied brain-PAD
DevAge_CB_BMI=lm(devAge_CB ~  Age + as.factor(Sex) + BMI_24 + as.factor(SCZ_PRS) + BMI_24*as.factor(SCZ_PRS), data = data)
summary(DevAge_CB_BMI)
ols_plot_resid_qq(DevAge_CB_BMI)
ols_test_normality(DevAge_CB_BMI)
ols_test_correlation(DevAge_CB_BMI)
ols_plot_resid_fit(DevAge_CB_BMI)
ols_plot_resid_hist(DevAge_CB_BMI)


# IV: PE by age 24 (binary)
# DV: ENGIMA-derived brain-PAD
DevAge_ENG_PE=lm(devAge_ENG ~ Age + Sex + as.factor(PE_24_18_bin) + as.factor(SCZ_PRS) + as.factor(SCZ_PRS)*as.factor(PE_24_18_bin), data = data)
summary(DevAge_ENG_PE)
# check linear regression assumptions
ols_plot_resid_qq(DevAge_ENG_PE)
ols_test_normality(DevAge_ENG_PE)
ols_test_correlation(DevAge_ENG_PE)
ols_plot_resid_fit(DevAge_ENG_PE)
ols_plot_resid_hist(DevAge_ENG_PE)

# DV: CentileBrain-derived brain-PAD
DevAge_CB_PE=lm(devAge_CB ~ Age + Sex + as.factor(PE_24_18_bin) + + as.factor(SCZ_PRS) + as.factor(SCZ_PRS)*as.factor(PE_24_18_bin), data = data)
summary(DevAge_CB_PE)
# check linear regression assumptions
ols_plot_resid_qq(DevAge_CB_PE)
ols_test_normality(DevAge_CB_PE)
ols_test_correlation(DevAge_CB_PE)
ols_plot_resid_fit(DevAge_CB_PE)
ols_plot_resid_hist(DevAge_CB_PE)

## done with exploratory analyses

## Done ##