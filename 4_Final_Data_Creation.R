#----------------- 4: Final Data Creation -----------------#

# This script is unfinished 

# This script performs the actions of Final Data Creation, defined by this flowchart: https://whimsical.com/tds-r-scripts-and-data-flow-VmAm6BzY1jUML2a32569t2
# It is designed to run on the HPC server

# It requires
# > PRS data
# > HES IBD Extraction
# > HES Other Outcomes Extraction
# > UKB Cancer Status
# > UB Covariates
# > Several packages
library(data.table)
library(openxlsx)
library(tidyverse)
library(tictoc)

args=commandArgs(trailingOnly=TRUE)
path=as.character(args[1])
print(path)
data<-paste0(path,"/data/")
setwd(data)
print("We've set the working directory")

print("This is the start of Step 4")
tic("Step 4")

# It outputs a single tidy dataset containing all the data we need for our analyses and validation

#----------------- 4a Mega-join of all the important columns -----------------#
# Original script is ukb_hes_join_everything.R by Xin

# Loading the data
print("Loading in the data now...")

# FROM UKB data
ukb_cancer_status<-readRDS("ukb_cancer_status.rds")
prs_all<-readRDS("PRS_all.rds")

# From UKB data
# Andrea's ukb_ML_covars.rds contains whats in Nas' earlier version
#nas_covar<-readRDS("ukb_covar.rds")
ukb_covar<-readRDS("ukb_ML_covars.rds")
ukb_drugs<-readRDS("ukb_drugs.rds")
ukb_biomarker<-readRDS("ukb_biomaker_final.rds")

# From HES data
# This covar does not seem to have removed withdrawn particpants
hes_covar<-readRDS("hes_covar.rds")
hes_ibd<-readRDS("hes_ibd_extraction.rds")
hes_other_outcomes<-readRDS("hes_other_outcomes.rds")

# Merge all UKB data ----------------------------------------------------
print("Now merging the data")

# We used a merge() here as PRS would be missing for some
# 10/03/2021 XZ: now left_join() for variable selection purposes
prs_all$eid<-as.integer(prs_all$eid)
ukb_cancer_prs<-left_join(ukb_cancer_status,prs_all,by="eid")
ukb_cancer_prs_covar<-merge(ukb_cancer_prs,ukb_covar,by="eid")
ukb_cancer_prs_covar_drugs<-merge(ukb_cancer_prs_covar,ukb_drugs,by="eid")
ukb_everything<-merge(ukb_cancer_prs_covar_drugs,ukb_biomarker,by="eid")

print(paste0("Number of participants in ukb cancer status: ", nrow(ukb_cancer_status)))
print(paste0("Number of participants in PRS score: ", nrow(prs_all)))
print(paste0("After left join: ", nrow(ukb_cancer_prs)))
print(paste0("After left join with covariates: ", nrow(ukb_cancer_prs_covar)))
print(paste0("After left join with medication: ", nrow(ukb_cancer_prs_covar_drugs)))
print(paste0("After left join with hes data: ", nrow(ukb_everything)))

# Merge all HES data ----------------------------------------------------

hes_ibd_covar<-merge(hes_ibd,hes_covar,by="eid")

hes_other_outcomes$eid<-row.names(hes_other_outcomes)
hes_everything<-merge(hes_ibd_covar,hes_other_outcomes,by="eid")


# Mega join -------------------------------------------------------------

ukb_hes_everything<-left_join(ukb_everything,hes_everything,by="eid") 

# Process diabetes from two different sources  --------------------------

ukb_hes_everything$diabtes_merge<-ifelse(ukb_hes_everything$Diabetes==1,1,
                                         ifelse(ukb_hes_everything$diabetes==1,1,0))

# Remove columns we don't need anymore (e.g. >33% missing data)
ukb_hes_everything <- ukb_hes_everything %>% 
  select(-c("birth_weight","duration_vigorous_activity",
            "nonbutter_spread_type","alcohol_with_meals","diabetes","Diabetes"))

print("Just saving the file ukb_hes_everything.rds ...")

# Export file 
saveRDS(ukb_hes_everything,"ukb_hes_everything.rds")

print("Step 4 is done. All data creation is now complete!!")

toc()