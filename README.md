# ASLPAC-SCZ-RbG-BrainAge analyses

This repository provides code/scripts for secondary data analyses as reported in the published manuscript: 

Constantinides, C., Baltramonaityte, V., Caramaschi, D., Han, L. K. M., Lancaster, T. M., Zammit, S., Freeman, T. P., & Walton, E. (2024). Assessing the association between global structural brain age and polygenic risk for schizophrenia in early adulthood: A recall-by-genotype study. Cortex; a journal devoted to the study of the nervous system and behavior, 172, 1–13. https://doi.org/10.1016/j.cortex.2023.11.015

Here we used data from the ALSPAC-Schizophrenia-Recall-by-Genotype (SCZ-RbG) imaging study as described previously (see data note by Sharp et al., 2020). 

### Step 1. Sample and variable selection

The following R script was used to extract the original ALSPAC-SCZ-RbG imaging sample (and variables of interest) from the wider ALSPAC birth cohort: 

(1) '1_Extract_ALSPAC_SCZ_RbG_Imaging_Dataset_20240928.R'

A desicription of the ALSPAC variables used for the current analyses can be found in the following Excel spreadsheet: 'Variables_Info_ALSPAC_SCZ_RbG_BA_CC_02102024.xlsx'. Note that column names for imaging-derived phenotype (IDPs) may differ in future data releases. See spreadsheet for column names as they currently appear in the ALSPAC data disctionary/variable search tool (https://variables.alspac.bris.ac.uk).

Imaging derived-phenotypes (IDPs) were previously derived from sMRI/T1-weighted scans (Sharp et al., 2020). 

### Step 2. Brain age prediction using the ENIGMA Brain Age model (Han et al., 2020) 
(NB. You may skip this step if you want to use the CentileBrain model only)

The following R scripts were used sequentially to prepare the model input files for the ALSPAC-SCZ-RbG imaging sample:

(2a) '2a_ENIGMA_Input_Prep_Step_2a_20240920.R' - formats column names based on the headers used by the ENIGMA model ('ENIGMA_headers.csv'). 

(2b) '2b_ENIGMA_Input_Prep_Step_2b_20240928.R' - takes the average of input features and splits dataset with respect to sex. 

Plus two additional scripts that are loaded by script 2b as functions:
'prepare.files.R',
'get.means.R'.

After running scripts 2a and 2b, the resulting model input files (males_raw.csv/females_raw.csv) should be uploaded on the PHOTONAI platform to get brain age predictions: 

https://photon-ai.com/enigma_brainage

Once you are on the platform, make sure that you click on the correct sex group (males/females) before uploading the corresponding input file. The two output files should be renamed to 'males_raw_out' and 'females_raw_out' for males and females respectively. 

Finally, the following script was run to generate model performance metrics/plots for the current sample and to prepare dataset for downstream statistical analyses (step 4):

(2c) '2c_ENIGMA_Output_Prep_Step_2c_20240928.R'

Plus two additional scripts that loaded by script 2c as functions: 'model.fits.brainAge.R', 
'model.fits.brainAge.SCZ_PRS.R'.

The ouput datasets from script 2c named as 'ALSPAC_SCZ_RbG_BA_ENG.RData' should be used for step 4a. 

### Step 3. Brain age prediction using the CentileBrain Brain Age model (Yu et al., 2024)

The following R script was used to prepare model input files for the ALSPAC-PE imaging sample:

(3a) '3a_CentileBrain_Input_prep_20240928.R' -  formats column names based on the headers used by CentileBrain model ('CentileBrain_headers.csv'). 

After running script 3a, the resulting model input files (CB_males_raw.xlsx/CB_females.raw.xlsx) should be uploaded on the CentileBrain platform to get brain age predictions: 

https://centilebrain.org/#/brainAge2

Once you are on the platform, make sure that you click on the correct age group (5≤age≤40 years) and sex group (males/females) before uploading the corresponding input file. You would just just need to download 'predicted age' for each sex group (top downlead button). The two output files should be renamed to 'CB_males_raw_out' and 'CB_females_raw_out' for males and females respectively. 

#### ** NOTE **
_At the time of accessing this newly developed resource for the current project (20.12.2022), the model used was trained on the age range of 20-30 years (for more details please refer to the supplementary material of Contantinides et al., 2024). As indicated above, this is NOT the same model as the one currrently available on the CentileBrain platofrom (link provided above), which was trained on the optimal age range of 5-40 years (refe to Yu et al., 2024 for the rationale and empirical data supporting this decision). Therefore, if you repeat the above process with the current version of the CentileBrain model for the current sample (ALSPAC-SCZ-RbG) you might get somewhat different results than those reported in the published manuscript (Constantinides et al., 2024)._
    
Finally, the following script was run to generate model performance metrics/plots for the current sample and to prepare the dataset for downstream statistical analyses (step 4):

(3b) '3b_CentileBrain_Output_prep_20240920.R'

Plus two additional scripts loaded by script 3b as functions: 'CB.model.fits.brainAge.R', 
'CB.model.fits.brainAge.SCZ_PRS.R'.

The ouput dataset named as 'ALSPAC_SCZ_RbG_BA_CB.RData' should be used for step 4. 

### Step 4. Statistical analyses

The following scripts were used for data wrangling and statistical analyses of the association between SCZ-PRS and brain-PAD: 

(4a) '4a_ALSPAC_SCZ_RbG_data_tranform_and_cleaning_28092024.R' - Recodes/transforms original variables and removes particpants with failed image QC

(4b) '4b_ALSPAC_SCZ_RbG_BA_data_analysis_28092024.R' - generates sample descriptives (including model pefromance for the final sample) and statistical ouput for assocation between SCZ-PRS and brain-PAD

-End of workflow-

### References

Constantinides, C., Baltramonaityte, V., Caramaschi, D., Han, L. K. M., Lancaster, T. M., Zammit, S., Freeman, T. P., & Walton, E. (2024). Assessing the association between global structural brain age and polygenic risk for schizophrenia in early adulthood: A recall-by-genotype study. Cortex; a journal devoted to the study of the nervous system and behavior, 172, 1–13. https://doi.org/10.1016/j.cortex.2023.11.015

Han, L. K. M., Dinga, R., Hahn, T., Ching, C. R. K., Eyler, L. T., Aftanas, L., Aghajani, M., Aleman, A., Baune, B. T., Berger, K., Brak, I., Filho, G. B., Carballedo, A., Connolly, C. G., Couvy-Duchesne, B., Cullen, K. R., Dannlowski, U., Davey, C. G., Dima, D., Duran, F. L. S., … Schmaal, L. (2021). Brain aging in major depressive disorder: results from the ENIGMA major depressive disorder working group. Molecular psychiatry, 26(9), 5124–5139. https://doi.org/10.1038/s41380-020-0754-0

Sharp, T. H., McBride, N. S., Howell, A. E., Evans, C. J., Jones, D. K., Perry, G., Dimitriadis, S. I., Lancaster, T. M., Zuccolo, L., Relton, C., Matthews, S. M., Breeze, T., David, A. S., Drakesmith, M., Linden, D. E. J., Paus, T., & Walton, E. (2020). Population neuroimaging: generation of a comprehensive data resource within the ALSPAC pregnancy and birth cohort. Wellcome open research, 5, 203. https://doi.org/10.12688%2Fwellcomeopenres.16060.1

Yu, Y., Cui, H. Q., Haas, S. S., New, F., Sanford, N., Yu, K., Zhan, D., Yang, G., Gao, J. H., Wei, D., Qiu, J., Banaj, N., Boomsma, D. I., Breier, A., Brodaty, H., Buckner, R. L., Buitelaar, J. K., Cannon, D. M., Caseras, X., Clark, V. P., … ENIGMA‐Lifespan Working Group (2024). Brain-age prediction: Systematic evaluation of site effects, and sample age range and size. Human brain mapping, 45(10), e26768. https://doi.org/10.1002/hbm.26768
