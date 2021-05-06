#----------------- 6.1: Sensitivity Analysis & Assumptions Check -----------------#

# This script is complete and takes about XXXX mins to run

# This script performs the actions of Sensivity Analysis & Assumptions Check, defined by this flowchart: https://whimsical.com/tds-r-scripts-and-data-flow-VmAm6BzY1jUML2a32569t2
# It is designed to run on the HPC server

# It requires
# > all previous steps
# > UKB genotype data
# > Several packages
library(tidyverse)#installed
library(tictoc)#installed
library(data.table)#installed
library(MendelianRandomization)
# Documentation: https://cran.r-project.org/web/packages/MendelianRandomization/MendelianRandomization.pdf
# Github repo: https://github.com/cran/MendelianRandomization



args=commandArgs(trailingOnly=TRUE)
nchunks=as.numeric(args[1])
path=as.character(args[2])
print(path)
datapath<-paste0(path,"/data/")
analysis<-paste0(path,"/analysis/")
setwd(datapath)
print("We've set the working directory")

print("Welcome to the start of the Job, we are initiating Step 6 now (Sensitivity Analysis)...")
tic("Step 6.1")

# It outputs analysis and visualisation in the folder "analysis"

#----------------- 6a Prepare fam file and List IBD -----------------#
# Original code is 0-2-make_fam_snps.R by Xin
# This step works and takes between XX & XX mins to run
print("Starting 6a")
tic("Step 6a")

# Read in libraries

# Set directory
dest<-datapath
output<-analysis

setwd(dest)

# Read in data for case status
# use the mega file for consistency with other ongoing workstreams
covar<-readRDS("ukb_hes_everything.rds")

# Remove participants where PRS_ibd is not available
covar<-covar[!is.na(covar$PRS_ibd),]

# Set path for imputed ukb genotype data
#ukb_imp="/rds/general/project/hda_students_data/live/Group4/General/"

##################

# Create a vector of cases
case=rep(1,nrow(covar))  # 1 for controls
case[which(covar$colon==1)]=2  # 2 for cases
names(case)=covar$eid
table(case) # rough check

fam=data.frame(fread("ukb_imp.fam"))
rownames(fam)=fam[,1]
ids=intersect(rownames(fam), names(case))
fam[,6]=0
fam[ids,6]=case[ids]
print(table(fam[,6]))

write.table(fam, "0-2-ukb_imp_yayy.fam", row.names=FALSE, col.names=FALSE, quote=FALSE)

##### Extract SNPs list

snps=read.csv("iv_maf_ibd.txt")$rsid

write.table(snps,paste0(dest,"0-2-ibd_snps_list.txt"), row.names=FALSE,col.names=FALSE, quote=FALSE)


toc()

#----------------- 6b Prepare confounders for mini GWAS IBD -----------------#
# Original code is 0-1-prep_covars.R by Xin
# This step works and takes between XX & XX mins to run
print("Starting 6b")
tic("Step 6b")

# Define directory
dest<-datapath
ukb<-"ukb26390.csv"

# Load and prepare the variables
# use the mega file for consistency with other ongoing workstrams
covar<-readRDS("ukb_hes_everything.rds")

# Remove participants where PRS_ibd is not available
covar<-covar[,c("eid","age","gender","PRS_ibd")]
covar<-covar[!is.na(covar$PRS_ibd),]
covar$gender<-as.factor(covar$gender)

# Selected eid as required by the original file
# Age and gender are confounders required ror SNP-level GWAS [checked with Sonja]
covar<-covar[,c("eid","age","gender")]

# Extract PCs provided by Barbara
pc<-readRDS("GWAS_PCs.rds")
pc$eid<-as.integer(pc$eid) # extract row name of PC for eid
pc<-pc[,1:11]  # keep first 10 PCS only, first column is eid

# Join confounders with 10 PCs
covar_gwas<-left_join(covar,pc,by="eid")

### Preparing variables ------------------------------------------------

# covar$eid=as.numeric(covar$eid)  # not required as already integer
# covar$age=as.numeric(covar$age)  # not required as already integer

#covars$immuno=as.factor(covars$immuno) # not required
#covars$Batch=as.factor(as.character(covars$Batch)) # not required
conf_list=colnames(covar_gwas)[]
#print(conf_list)

mymodel=model.matrix(as.formula(paste0("~", paste(conf_list, collapse="+"))), data=covar_gwas)
#print(head(mymodel))
mymodel=mymodel[,-1]
colnames(mymodel)[1]="IID"
mymodel=cbind(FID=mymodel[,1], mymodel)
colnames(mymodel)=gsub(" ", "", colnames(mymodel))
#print(mymodel)

# Saving confounders file
write.table(mymodel, "0-1-confounders.txt", row.names=FALSE, col.names=TRUE, quote=FALSE, sep=" ")

toc()

#----------------- 6c Extract Phenoscanner output for IBD SNPs -----------------#
# Original code is 1-2-phenoscanner.R by Xin
# This step works and takes between XX & XX mins to run
print("Starting 6c")
tic("Step 6c")

# Define paths and set working directory
y_path<-datapath
X_path<-datapath

# Set working directory

### Data prep -----------------------------------------------------

# Read in data
X<-readRDS("corrected_betas_ibd.rds")  # beta=logOR and corrected, no se

# Extract list of SNPs for later
X_snps<-X$id

# Pheno scanner -----------------------------------------------------
print("Starting phenoscanner section")
# Export phenoscanner output to file
for (i in (1:length(X_snps))){  # XZ updated 05/05/2021
  
  snp<-X_snps[i] # XZ updated 05/05/2021
  snp_data<-phenoscanner(snp)
  write.table(snp_data,"1-2-phenoscanner_output.csv",
              append=TRUE,sep="/")
  
  print(paste0("Progress: ",i,"/",length(X_snps)))  
}
toc()
dir.create(file.path(getwd(), "v4"), showWarnings = FALSE)

#----------------- 6d Run mini GWAS -----------------#
# Original code is 0-3-logistic_geno.sh by Xin
##### This is done in an sh script!

print("Step 6.1 complete! Now for the mini GWAS.")
toc()