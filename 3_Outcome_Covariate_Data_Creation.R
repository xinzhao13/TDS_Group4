#----------------- 3: Outcome and Covariate Data Creation -----------------#

# This script performs the actions of Outcome and Covariate Data Creation, defined by this flowchart: https://whimsical.com/tds-r-scripts-and-data-flow-VmAm6BzY1jUML2a32569t2
# It is designed to run on the HPC server

# It requires
# > HES Data
# > Participant withdrawal data
# > UKB Participant Data (ukb2390.csv)
# > UKB Biomarker Data (ukb27725.csv)
# > ukb_drugs_list.txt by Marie
# > Several packages
library(data.table)
library(openxlsx)
library(tidyverse)
library(bit64)
library(dplyr)
library(tictoc)

args=commandArgs(trailingOnly=TRUE)
path=as.character(args[1])
print(path)
data<-paste0(path,"/data/")
setwd(data)
print("We've set the working directory")

path_data<-data
path_out<-data
path_list<-data

print("This is the start of Step 3")
tic("Step 3")

# It outputs seven files with outcome and covariate data for participants

#----------------- 3a Get Biomarker Data -----------------#
# Original script is biomarker_extraction_final.R by Nas
print("Step 3a starting now")
tic("Step 3a")

# Loading the data
mydata=data.frame(fread("ukb27725.csv", nrows=5))
myfields=unname(unlist(read.table("List_field_ids_to_extract_biomarker.txt", header=FALSE)))

# Read withdrawn eid data----------
withdrawn_eid<-read.csv("w19266_20200204.csv")[,1]

ukb_eid<-fread("ukb26390.csv", select=1)$eid

ukb_consent_eid<-ukb_eid[!(ukb_eid %in% withdrawn_eid)]

## -------------------------------------------

# Extracting the column ids 
column_id=grep("eid", colnames(mydata))
found_fieldids=NULL
for (k in 1:length(myfields)){
  mygrep=grep(paste0("X",myfields[k],"."), fixed=TRUE, colnames(mydata))
  if (length(mygrep)>0){
    found_fieldids=c(found_fieldids, myfields[k])
  }
  
  column_id=c(column_id, mygrep)
}
#View(myfields)
# Extracting required columns from dataset
extracted=data.frame(fread("ukb27725.csv", select=column_id))

# Remove withdrawn participants
extracted<-filter(extracted,eid %in% ukb_consent_eid)

## -------------------------------------------

biomarker <- extracted
#View(biomarker)

# Assign value to basline if there is NAs
for(r in 1:nrow(biomarker) ) {
  if(is.na(biomarker[r,2]) && !is.na(biomarker[r,3]) ) {
    for(c in 2:(ncol(biomarker)-1) ) {
      #print(cbind(r,c) )
      if(c%%2 == 0) {
        biomarker[r,c] <- biomarker[r,c+1]
      }
    }
  }
}

#vis_miss(biomarker[1:100,])

# Get column names
#colnames(biomarker)

#Remove excess columns
biomarker <- select(biomarker, -ends_with("1.0"))
#View(biomarker)

# Remove Oestradiol & Rheumatoid factor
drop <- c("X30800.0.0","X30820.0.0")
biomarker_clean = biomarker[,!(names(biomarker) %in% drop)]

#nrow(biomarker)

#Remove NAs
#biomarker2 <- biomarker_clean[complete.cases(biomarker_clean), ]
#nrow(biomarker2)

# Modify column names
colnames(biomarker_clean) <- c("eid","alanine_aminotransferase", "albumin", "alkaline_phosphatase", "apolipoprotein_A", "apolipoprotein_B", "aspartate_aminotransferase", "C-reactive_protein", "calcium", "cholesterol", "creatinine", "cystatin_C", "direct_bilirubin", "gamma_glutamyltransferase", "glucose", "glycated_haemoglobin(HbA1c)", "hdl_cholesterol", "igf-1", "ldl_direct", "lipoprotein_A", "phosphate", "shbg", "testosterone", "total_bilirubin", "total_protein", "triglycerides", "urate", "urea", "vitamin_D" ) 

#View(biomarker2)
#summary(biomarker2)
#View(ukbiomaker)

# Save RDS for ukbiomarker_final.rds
saveRDS(biomarker_clean, "ukb_biomaker_final.rds")

toc()

#----------------- 3b Get HES Covariate data -----------------#
# Original script is hes_covar_extraction.R by Xuyi
print("Beginning Step 3b")
tic("Step 3b")

### Creating the dataset--------------

mydata=data.frame(fread("ukb26390.csv", select=1))
rownames(mydata)=mydata[,1] #rowname is eid

diseases=c("Diabetes","unspecific_benign_neoplasm_of_colon","FH_MNDO")
#unspecific_benign_neoplasm_of_colon:Adenomatosis of colon,Large intestine NOS,Polyposis (hereditary) of colon
#FH_MNDO:Family history of malignant neoplasm of digestive organs

names(diseases)=diseases
x=matrix(0,nrow=nrow(mydata),ncol=length(diseases))
colnames(x)=names(diseases)
mydata=cbind(mydata,x)
hes=data.frame(fread("hesin_diag.txt"))

### Disease definition------------------

# <--- below are all user-input
# naming convention is icd10_[disease_name]

icd10_Diabetes=c(paste0("E", 10:14)) 
icd10_unspecific_benign_neoplasm_of_colon="D126"
icd10_FH_MNDO = "Z800"

### Build icd_list for diseases of interest--------

for (d in 1:length(diseases)){
  print(names(diseases)[d])
  disease=diseases[d]
  icd10_list=eval(parse(text=paste0("icd10_",disease)))
  myeids=NULL
  
  # ICD10 codes
  pb=txtProgressBar(style=3)
  for (k in 1:length(icd10_list)){
    setTxtProgressBar(pb, k/length(icd10_list))
    tmp=as.character(hes$eid[grepl(paste0("^", icd10_list[k]), hes$diag_icd10)])
    myeids=c(myeids, tmp)
  }
  myeids=unique(myeids)
  cat("\n")
  print(length(myeids))
  
  print(table(mydata[myeids,names(diseases)[d]]))
  mydata[myeids,names(diseases)[d]]=1
  print(table(mydata[,names(diseases)[d]]))
  cat("\n")
}

# Save as RDS file
saveRDS(mydata, "hes_covar.rds")

toc()


#----------------- 3c Get HES non-cancer outcomes -----------------#
# Original script is hes_otheroutcome_extraction.R by Xuyi
print("Beginning Step 3c")
tic("Step 3c")

### Creating the dataset--------------

mydata=data.frame(fread("ukb26390.csv", select=1))
rownames(mydata)=mydata[,1] #rowname is eid

diseases=c("Depressive_episode","Recurrent_depressive_disorder",
           "Phobic_anxiety_disorders","Other_anxiety_disorders")

names(diseases)=diseases
x=matrix(0,nrow=nrow(mydata),ncol=length(diseases))
colnames(x)=names(diseases)
mydata=cbind(mydata,x)
hes=data.frame(fread("hesin_diag.txt"))

### Disease definition------------------

# <--- below are all user-input
# naming convention is icd10_[disease_name]

icd10_Depressive_episode="F32" 
icd10_Recurrent_depressive_disorder="F33" 

icd10_Phobic_anxiety_disorders="F40"
icd10_Other_anxiety_disorders="F41"

### Build icd_list for diseases of interest--------

for (d in 1:length(diseases)){
  print(names(diseases)[d])
  disease=diseases[d]
  icd10_list=eval(parse(text=paste0("icd10_",disease)))
  myeids=NULL
  
  # ICD10 codes
  pb=txtProgressBar(style=3)
  for (k in 1:length(icd10_list)){
    setTxtProgressBar(pb, k/length(icd10_list))
    tmp=as.character(hes$eid[grepl(paste0("^", icd10_list[k]), hes$diag_icd10)])
    myeids=c(myeids, tmp)
  }
  myeids=unique(myeids)
  cat("\n")
  print(length(myeids))
  
  print(table(mydata[myeids,names(diseases)[d]]))
  mydata[myeids,names(diseases)[d]]=1
  print(table(mydata[,names(diseases)[d]]))
  cat("\n")
}

# Combine Depressive_episode and Recurrent_depressive_disorder into depression
mydata_other<-mydata %>% 
  mutate(depression=ifelse(Depressive_episode == 1 | Recurrent_depressive_disorder == 1, 1, 0),
         anxiety=ifelse(Phobic_anxiety_disorders == 1 | Other_anxiety_disorders == 1, 1, 0))

# Save as RDS file
saveRDS(mydata_other, "hes_other_outcomes.rds")

toc()


#----------------- 3d Get HES data for Validation -----------------#
# Original script is hes_ibd_extraction.R by Xuyi
print("Now we're onto Step 3d")
tic("Step 3d")

### Creating the dataset--------------

mydata=data.frame(fread("ukb26390.csv", select=1))
rownames(mydata)=mydata[,1] #rowname is eid

diseases=c("cd","uc")

names(diseases)=diseases
x=matrix(0,nrow=nrow(mydata),ncol=length(diseases))
colnames(x)=names(diseases)
mydata=cbind(mydata,x)
hes=data.frame(fread("hesin_diag.txt"))

### Disease definition------------------

# <--- below are all user-input
# naming convention is icd10_[disease_name]

#icd10_crc="C18" # colon cancer 
icd10_cd="K50" # crohn's disease
icd10_uc="K51" # ulcerative colitis

### Build icd_list for diseases of interest--------

for (d in 1:length(diseases)){
  print(names(diseases)[d])
  disease=diseases[d]
  icd10_list=eval(parse(text=paste0("icd10_",disease)))
  myeids=NULL
  
  # ICD10 codes
  pb=txtProgressBar(style=3)
  for (k in 1:length(icd10_list)){
    setTxtProgressBar(pb, k/length(icd10_list))
    tmp=as.character(hes$eid[grepl(paste0("^", icd10_list[k]), hes$diag_icd10)])
    myeids=c(myeids, tmp)
  }
  myeids=unique(myeids)
  cat("\n")
  print(length(myeids))
  
  print(table(mydata[myeids,names(diseases)[d]]))
  mydata[myeids,names(diseases)[d]]=1
  print(table(mydata[,names(diseases)[d]]))
  cat("\n")
}

# Combine uc and cd into ibd
mydata_ibd<-mydata %>% 
  mutate(ibd=ifelse(uc == 1 | cd == 1, 1, 0))

# Save as RDS file
saveRDS(mydata_ibd, "hes_ibd_extraction.rds")

toc()

#----------------- 3e Get UKB Cancer Outcomes for All Chosen Cancers -----------------#
# Original script is ukb_cancer_extraction.R by Xin
print("Beginning Step 3e")
tic("Step 3e")

### Define path ------------------------------------------------------------------

# Loading the data

### Trim UKB with only 40006 columns + remove withdrawn participants -------------

# Read eid data
withdrawn_eid<-read.csv("w19266_20200204.csv")[,1]
ukb_eid<-fread("ukb26390.csv", select=1)$eid

ukb_consent_eid<-ukb_eid[!(ukb_eid %in% withdrawn_eid)]

# Read ukb header
ukb<-data.frame(fread("ukb26390.csv",nrows=0))

# Check for ICD10 data field
colnames(ukb)[grep("40006",colnames(ukb))] # 40006=cancer ICD10

# Get column index
colids<-c(1, grep("40006",colnames(ukb)))

# Read the actual ukb data + selected columns
ukb_cancer<-data.frame(fread("ukb26390.csv",select=colids))

# Remove withdrawn participants
ukb_cancer_consent<-filter(ukb_cancer,ukb_eid %in% ukb_consent_eid)

# Save cancer status of all UKB participants with consent
#saveRDS(ukb_cancer_consent,paste0(path_out,"ukb_cancer_consent.rds"))

# Print message
print(paste0(length(ukb_eid)-length(ukb_consent_eid),
             " participants removed"))

### Extract cancer status -------------------------------------------

# Read cancer status of all UKB participants with consent
#ukb_cancer_consent<-readRDS(paste0(path_out,"ukb_cancer_consent.rds"))
ukb_consent_eid<-ukb_cancer_consent$eid

# Define cancer types of interest
cancers<-c("C17", # small intestine
           "C18", # all colon cancer types
           "C19", # rectosigmoid junction
           "C20", # rectal cancer
           "C21", # anal
           "C22") # liver

# Define names
cancer_names<-c("small_intestine","colon","rectosigmoid","rectal","anal","liver")

# Create empty df
ukb_cancer_status_all<-data.frame(matrix(NA,nrow=nrow(ukb_cancer_consent),ncol=1))
colnames(ukb_cancer_status_all)<-c("eid")
ukb_cancer_status_all$eid<-ukb_cancer_consent$eid

# Create cancer status by cancer type
for (i in 1:length(cancers)){
  
  cancer<-cancers[i]
  cancer_name<-cancer_names[i]
  
  print(paste0("Start: ",cancer))
  
  # Filter C180-C189 columns to identify participants with colon cancer
  ukb_cancer_cases<-filter_all(ukb_cancer_consent, any_vars(str_detect(., cancer)))
  
  # Extract list of eids by cc case and control
  ukb_cancer_cases_eid<-ukb_cancer_cases$eid
  ukb_cancer_controls_eid<-ukb_cancer_consent$eid[!(ukb_cancer_consent$eid %in% ukb_cancer_cases_eid)]
  
  # Check if match
  print("Check: number of cancer cases + controls = number of consent")
  print(length(ukb_cancer_cases_eid)+length(ukb_cancer_controls_eid)==length(ukb_consent_eid))
  
  # Combine into dataframe with eid and cc status (1=case, 0=control)
  ukb_cancer_status<-data.frame(matrix(0,nrow=nrow(ukb_cancer_consent),ncol=2))
  colnames(ukb_cancer_status)<-c("eid","cancer_status")
  
  ukb_cancer_status$eid<-ukb_cancer_consent$eid
  ukb_cancer_status$cancer_status<-ifelse(ukb_cancer_status$eid %in% ukb_cancer_cases_eid,1,0)
  
  colnames(ukb_cancer_status)<-c("eid",paste0(cancer_name))
  
  ukb_cancer_status_all<-merge(ukb_cancer_status_all,ukb_cancer_status,by="eid")
  
}

ukb_cc_status<-ukb_cancer_status_all[,c("eid","colon")]
colnames(ukb_cc_status)<-c("eid","cc_status")

# Export into rds file
#saveRDS(ukb_cc_status, paste0(path_out,"ukb_cc_status.rds"))
saveRDS(ukb_cancer_status_all, paste0(path_out,"ukb_cancer_status.rds"))

toc()

#----------------- 3f Get Covariates for ML Confounder Analysis -----------------#
# Original script is ukb_ML_covars.R by Andy
print("Starting Step 3f")
tic("Step 3f")

##-----------barbara code to get the data we want efficiently
# Loading the data
mydata=data.frame(fread("ukb26390.csv", nrows=5))
myfields=unname(unlist(read.table("List_field_ids_to_extract.txt", header=FALSE)))

# Extracting the column ids 
column_id=grep("eid", colnames(mydata))
found_fieldids=NULL
for (k in 1:length(myfields)){
  mygrep=grep(paste0("X",myfields[k],"."), fixed=TRUE, colnames(mydata))
  if (length(mygrep)>0){
    found_fieldids=c(found_fieldids, myfields[k])
  }
  column_id=c(column_id, mygrep)
}

# Extracting required columns from dataset
extracted=data.frame(fread("ukb26390.csv", select=column_id)) # Path to change!
# saveRDS(extracted, "/rds/general/project/hda_students_data/live/Group4/General/andrea/ML_covars_preadjustment_tobedeleted.rds")
##-----------end of barbara code
print("Barbara Code Finished")

# extracted<- readRDS("/rds/general/project/hda_students_data/live/Group4/General/andrea/ML_covars_preadjustment_tobedeleted.rds")

# Import the field IDs I want to use and convert them to the correct column format
## myfields=unname(unlist(read.table("/rds/general/project/hda_students_data/live/Group4/General/andrea/List_field_ids_to_extract.txt", header=FALSE)))
## myfields=as.character(myfields)
## preparefunction <- function(x) {
#1 Add an X to the start
#2 Add .0.0, to the end
##  paste("X",x,".0.0,", sep="", collapse=NULL)
## }
## myfields <- sapply(myfields, preparefunction)
## selectcols <- paste(myfields, collapse=" ")
# Manually copying the code from the above into the below to save time

# Select the columns I want to use
# extracted <- extracted %>% select(eid, X21022.0.0, X31.0.0, X21001.0.0, X20116.0.0, X738.0.0, X21000.0.0, X20022.0.0, X20119.0.0, X2443.0.0, X894.0.0, X914.0.0, X874.0.0, X981.0.0, X943.0.0, X991.0.0, X971.0.0, X884.0.0, X904.0.0, X864.0.0, X924.0.0, X1110.0.0, X1120.0.0, X10749.0.0, X10016.0.0, X1130.0.0, X1140.0.0, X1160.0.0, X1170.0.0, X1180.0.0, X1190.0.0, X1200.0.0, X20160.0.0, X10895.0.0, X1239.0.0, X1249.0.0, X2644.0.0, X3436.0.0, X3446.0.0, X5959.0.0, X3456.0.0, X6194.0.0, X6183.0.0, X3466.0.0, X3476.0.0, X3486.0.0, X3496.0.0, X3506.0.0, X6158.0.0, X2867.0.0, X2877.0.0, X2887.0.0, X2897.0.0, X2907.0.0, X10827.0.0, X6157.0.0, X10115.0.0, X2926.0.0, X2936.0.0, X1259.0.0, X1269.0.0, X1279.0.0, X1289.0.0, X1299.0.0, X1309.0.0, X1319.0.0, X1329.0.0, X1339.0.0, X1349.0.0, X1359.0.0, X1369.0.0, X1379.0.0, X1389.0.0, X3680.0.0, X6144.0.0, X10855.0.0, X1408.0.0, X1418.0.0, X1428.0.0, X2654.0.0, X10767.0.0, X1438.0.0, X1448.0.0, X10776.0.0, X1458.0.0, X1468.0.0, X1478.0.0, X1488.0.0, X1498.0.0, X1508.0.0, X1518.0.0, X1528.0.0, X1538.0.0, X1548.0.0, X10912.0.0, X20117.0.0, X1558.0.0, X3731.0.0, X4407.0.0, X4418.0.0, X4429.0.0, X4440.0.0, X4451.0.0, X4462.0.0, X1568.0.0, X1578.0.0, X1588.0.0, X1598.0.0, X1608.0.0, X5364.0.0, X1618.0.0, X1628.0.0, X2664.0.0, X10818.0.0, X3859.0.0, X10853.0.0, X2129.0.0, X2139.0.0, X2149.0.0, X2159.0.0, X3669.0.0)

# Rename all columns
extracted <- extracted %>% 
  transmute(eid,
            age = X21022.0.0, 
            gender = X31.0.0, 
            bmi = X21001.0.0, 
            smoking_status = X20116.0.0, 
            household_income = X738.0.0, 
            ethnicity = X21000.0.0, 
            birth_weight = X20022.0.0, 
            #employment_status = X20119.0.0, 
            diabetes = X2443.0.0, 
            
            # Physical activity
            duration_moderate_activity = X894.0.0, 
            duration_vigorous_activity = X914.0.0, 
            duration_walks = X874.0.0, 
            duration_pleasure_walking = X981.0.0, 
            freq_stair_climbing = X943.0.0, 
            #freq_strenuous_sports = X991.0.0, 
            freq_pleasure_walk = X971.0.0, 
            days_moderate_activity_per_week = X884.0.0, 
            days_vigorous_activity_per_week = X904.0.0, 
            days_walk_per_week = X864.0.0, 
            usual_walking_pace = X924.0.0, 
            
            # Electronics
            length_phone_use = X1110.0.0, 
            weekly_phone_use = X1120.0.0, 
            #time_using_phone = X10749.0.0, 
            #regular_use_handsfree_device = X10016.0.0, 
            handsfree_device_last3mo = X1130.0.0, 
            diff_in_phone_use = X1140.0.0, 
            
            # Sleep
            sleep_duration = X1160.0.0, 
            get_up_morning = X1170.0.0, 
            chronotype = X1180.0.0, 
            nap_during_day = X1190.0.0, 
            sleeplessness = X1200.0.0, 
            
            # Smoking
            ever_smoked = X20160.0.0, 
            #light_smokers_pilot = X10895.0.0, 
            current_tobacco_smoking = X1239.0.0, 
            past_tobacco_smoking = X1249.0.0, 
            #light_smokers_in_lifetime = X2644.0.0, 
            #age_started_smoking = X3436.0.0, 
            #type_tobacco_smoked = X3446.0.0, 
            #previously_smoked_most_days = X5959.0.0, 
            #number_cigarettes_currently = X3456.0.0, 
            #age_stopped_smoking_cigarettes = X6194.0.0, 
            #number_cigarettes_previously_smoked= X6183.0.0, 
            #time_waking_to_cigarette = X3466.0.0, 
            #difficulty_not_smoking_for_one_day = X3476.0.0, 
            #ever_tried_to_stop_smoking = X3486.0.0, 
            #want_to_stop_smoking = X3496.0.0, 
            #smoking_compared_to_10yrs = X3506.0.0, 
            #why_reduced_smoking = X6158.0.0, 
            #age_started_smoking_former_smokers = X2867.0.0, 
            #type_tobacco_previously_smoked = X2877.0.0, 
            #number_cigarettes_previously = X2887.0.0, 
            #age_stopped_smoking = X2897.0.0, 
            #ever_stopped_sixmo = X2907.0.0, 
            #ever_stopped_sixmo_pilot = X10827.0.0, 
            #why_stopped_smoking = X6157.0.0, 
            #why_stopped_smoking_pilot = X10115.0.0, 
            #number_unsuccessful_stopsmoking_attempts = X2926.0.0, 
            #likelihood_resume = X2936.0.0, 
            smokers_in_house = X1259.0.0, 
            exposure_tobacco_smoke_home = X1269.0.0, 
            exposure_tobacco_smoke_outside = X1279.0.0, 
            
            # Diet
            cooked_veg = X1289.0.0, 
            salad_raw_veg = X1299.0.0, 
            fresh_fruit = X1309.0.0, 
            dried_fruit = X1319.0.0, 
            oily_fish = X1329.0.0, 
            nonoily_fish = X1339.0.0, 
            processed_meat = X1349.0.0, 
            poultry = X1359.0.0, 
            beef = X1369.0.0, 
            lamb = X1379.0.0, 
            pork = X1389.0.0, 
            #age_last_ate_meat = X3680.0.0, 
            never_eat_egg_dairy_wheat_sugar = X6144.0.0, 
            #never_eat_pilot = X10855.0.0, 
            cheese = X1408.0.0,
            milk_type = X1418.0.0, 
            spread_type = X1428.0.0, 
            nonbutter_spread_type = X2654.0.0, 
            #spread_type_pilot = X10767.0.0, 
            bread_intake = X1438.0.0, 
            bread_type = X1448.0.0, 
            #bread_typeintake_pilot = X10776.0.0, 
            cereal = X1458.0.0, 
            cereal_type = X1468.0.0, 
            salt_added = X1478.0.0, 
            tea = X1488.0.0, 
            coffee = X1498.0.0, 
            coffee_type = X1508.0.0, 
            hot_drink_temp = X1518.0.0, 
            water = X1528.0.0, 
            major_diet_changes = X1538.0.0, 
            variation_in_diet = X1548.0.0, 
            #variation_in_diet_pilot = X10912.0.0,
            
            # Alcohol
            alcohol_drinker_status = X20117.0.0, 
            alcohol_intake_freq = X1558.0.0, 
            #former_alcohol_drinker = X3731.0.0, 
            #redwine_monthly = X4407.0.0, 
            #champagnewhitewine_monthly = X4418.0.0, 
            #beercider_monthly = X4429.0.0, 
            #spirits_monthly = X4440.0.0, 
            #fortifiedwine_monthly = X4451.0.0, 
            #otheralcohol_monthly = X4462.0.0, 
            redwine_weekly = X1568.0.0, 
            champagnewhitewine_weekly = X1578.0.0, 
            beercide_weekly = X1588.0.0, 
            spirits_weekly = X1598.0.0, 
            fortifiedwine_weekly = X1608.0.0, 
            #otheralcohol_weekly = X5364.0.0, 
            alcohol_with_meals = X1618.0.0, 
            alcohol_vs_10yrs = X1628.0.0, 
            #reason_for_reduction = X2664.0.0, 
            #reason_for_reduction_pilot = X10818.0.0, 
            #reason_stopped_drinking = X3859.0.0, 
            #reason_stopped_drinking_pilot = X10853.0.0, 
            
            # Sexual factors
            answered_sexual_history_qs = X2129.0.0, 
            age_secual_intercourse = X2139.0.0, 
            lifetime_sexual_partners = X2149.0.0, 
            ever_same_sex = X2159.0.0)
#lifetime_same_sex_sexual_partners = X3669.0.0)

print("Columns renamed")

# Read withdrawn eid data------------
withdrawn_eid<-read.csv("w19266_20200204.csv")[,1]

ukb_eid<-fread("ukb26390.csv", select=1)$eid

ukb_consent_eid<-ukb_eid[!(ukb_eid %in% withdrawn_eid)]

# Remove withdrawn participants
extracted<-filter(extracted,eid %in% ukb_consent_eid)

# Save file before we actually go ahead and adjust the variables to all be meaningful
# saveRDS(extracted, "/rds/general/project/hda_students_data/live/Group4/General/andrea/ML_covars_preadjustment_tobedeleted.rds")

print("Removed withdrawn participants")

print("Now we are beginning removing negatives with NAs")
## -------------------------------------------

# Now we need to check all the variables are coded correctly
# Change any miscoded values to NA
# We will not remove NAs - all data is retained at this stage (even though there are some very high NA counts!)

# -1, -3, -7, -10 and any negative number is indicative of an NA value (e.g. Do Not Know or Prefer Not to Answer)

# extracted <- readRDS("/rds/general/project/hda_students_data/live/Group4/General/andrea/ML_covars_preadjustment_tobedeleted.rds")
# summary(extracted)

colnames <- colnames(extracted)
for(i in 1:length(colnames)) {
  # extracted[,i] <- extracted[,i] %>% replace_with_na_all(condition = ~.x < 0)
  extracted[,i] <- ifelse( (!is.na(extracted[,i])) & (extracted[,i]<0) , NA , extracted[,i] )
  cat(colnames[i]," is DONE.\n")
}

# Some people have indicated 2000+ life sexual partners
# sum(ifelse(!is.na(lifetime_sexual_partners)&lifetime_sexual_partners>2000,1,0))
# I'm going to set these to NA
extracted$lifetime_sexual_partners <- ifelse( (!is.na(extracted$lifetime_sexual_partners)) & (extracted$lifetime_sexual_partners>2000) , NA , extracted$lifetime_sexual_partners )

# summary(extracted)

print("Replaced negative values with NAs")

print("Just need to save the file now")

#---------------------------------------------

# Save as RDS file
saveRDS(extracted, "ukb_ML_covars.rds")

toc()

#----------------- 3g Get drug covariates -----------------#
# Original script is ukb_drugs_extraction.R by Xin
print("Starting Step 3g")
tic("Step 3g")

### Define path ------------------------------------------------------------------

### Trim UKB to only 20003 columns -----------------------

# Set path
setwd(path_out)

# Read ukb header
ukb<-data.frame(fread(paste0(path_data,"ukb26390.csv"),nrows=0))

# Check for medication data field
colnames(ukb)[grep("20003",colnames(ukb))] # 20003=medication

# Get column index
colids<-c(1, grep("20003",colnames(ukb)))

# Read the actual ukb data + selected columns
ukb_drugs<-data.frame(fread(paste0(path_data,"ukb26390.csv"),select=colids))

### Remove withdrawn participants -----------------------

# Read eid data
withdrawn_eid<-read.csv(paste0(path_data,"w19266_20200204.csv"))[,1]
ukb_eid<-fread(paste0(path_data,"ukb26390.csv"), select=1)$eid

# Keep only particiapnts with consent
ukb_consent_eid<-ukb_eid[!(ukb_eid %in% withdrawn_eid)]

# Remove withdrawn participants
ukb_drugs_consent<-filter(ukb_drugs,ukb_eid %in% ukb_consent_eid)

# Convert/Ensure all columns except the first to/are characters
ukb_drugs_consent<-ukb_drugs_consent %>% mutate_at(-1,as.character)

# Save cancer status of all UKB participants with consent
#saveRDS(ukb_drugs_consent,paste0(path_out,"ukb_drugs_consent.rds"))

# Print message
print(paste0(length(ukb_eid)-length(ukb_consent_eid),
             " participants removed"))

### Extract medication data -------------------------------------------

# Read medication data of all UKB participants with consent
#ukb_drugs_consent<-readRDS(paste0(path_out,"ukb_drugs_consent.rds"))
ukb_consent_eid<-ukb_drugs_consent$eid

# Load list of drugs prepared by Marie
ukb_drugs_list<-read.csv(paste0(path_list,"ukb_drugs_list.txt"),
                         stringsAsFactors=FALSE)  # <-- to prevent R read data as factors

# Convert drug codes to characters (they were read as integers)
ukb_drugs_list$drug_code<-as.character(ukb_drugs_list$drug_code)

# Extract unique drug classes
drug_classes<-unique(ukb_drugs_list$drug_class)

# Create empty df
ukb_drugs_all<-data.frame(matrix(NA,nrow=nrow(ukb_drugs_consent),ncol=1))
colnames(ukb_drugs_all)<-c("eid")
ukb_drugs_all$eid<-ukb_drugs_consent$eid

# Create cancer status by cancer type
for (i in 1:length(drug_classes)){
  
  # Extract unique drug classes
  drug_class<-drug_classes[i]
  
  # Print drug class being processed
  print(paste0("Start: ",drug_class))
  
  # Extract drug names under class
  drug_codes<-ukb_drugs_list$drug_code[ukb_drugs_list$drug_class==drug_class]
  
  # Filter C180-C189 columns to identify participants with colon cancer
  ukb_drugs_sub<-filter_all(ukb_drugs_consent, any_vars(str_detect(., drug_codes)))
  
  # Extract list of eids taking current medication class
  ukb_drugs_yes_eid<-ukb_drugs_sub$eid
  
  # Combine into dataframe with eid and cc status (1=case, 0=control)
  ukb_drugs_sub<-data.frame(matrix(0,nrow=nrow(ukb_drugs_consent),ncol=2))
  colnames(ukb_drugs_sub)<-c("eid","drug_status")
  
  # Save medication status to df
  ukb_drugs_sub$eid<-ukb_drugs_consent$eid
  ukb_drugs_sub$drug_status<-ifelse(ukb_drugs_sub$eid %in% ukb_drugs_yes_eid,1,0)
  
  # Rename columns
  colnames(ukb_drugs_sub)<-c("eid",paste0(drug_class))
  
  # Inner join with main df
  ukb_drugs_all<-merge(ukb_drugs_all,ukb_drugs_sub,by="eid")
  
}

# Export into rds file
saveRDS(ukb_drugs_all, paste0(path_out,"ukb_drugs.rds"))
toc()

print("Step 3, Outcome & Covariate Data Creation, is now Complete!")
toc()
