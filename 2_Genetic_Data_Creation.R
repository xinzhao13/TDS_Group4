#----------------- 2: Genetic Data Creation -----------------#

# This script performs the actions of Genetic Data Creation
# It is designed to run on the HPC server

# It requires
# > IBD SNP and MAF Lists
# > UKB Genotyping Data ".bim" and ".fam" files
# > Several packages
library(snpStats)
library(genio)
library(parallel)
library(tictoc)

args=commandArgs(trailingOnly=TRUE)
nchunks=as.numeric(args[1])
path=as.character(args[2])
print(path)
data<-paste0(path,"/data/")
setwd(data)
print("We've set the working directory")

print("Starting Step 2 (Genetic Data Creation) R Script now...")
print("Step 2a has already been completed.")
tic("Step 2")

# It outputs a Polygenic Risk Score (PRS) of IBD, CD, and UC for all participants

#----------------- 2a Import SNPs -----------------#
# Here I need to integrate the sh scripts into my central sh script

# import_snps.sh
# merge_genetic_data_all.sh

# This will output a lovely file called UKB_imp_merged.bim within full_scripts/data/import_snps_results_all

#----------------- 2b UKB Genetic Data Creation -----------------#
# Original code is 3-read_from_plink.R
# This takes about 90 seconds
print("Starting Step 2b")
tic("Step 2b")

# Define path
res<-"import_snps_results_all/"

# Extract SNP names
bim=read.table(paste0(res,"ukb_imp_merged.bim"), stringsAsFactors=FALSE)
mysnps=bim[!duplicated(bim[,2]),2]

# Rename participants
fam=read.table("ukb_imp.fam")
fam[,1]=fam[,2]=1:nrow(fam)
write.table(fam, paste0(res,"ukb_imp_merged.fam"), quote=FALSE, row.names=FALSE, col.names=FALSE)

print("Reading some plink data now")
# Read plink data
mydata=read.plink(paste0(res,"ukb_imp_merged"), select.snps=mysnps, na.strings="-9")

# Extract genotype data
genodata=mydata$genotypes
genodata=data.frame(genodata)

print("Renaming participants with silly names")
# Rename participants
fam=read.table("ukb_imp.fam")
rownames(genodata)=fam[,1]

# Save genetic data
saveRDS(genodata, paste0(res,"genetic_data_extracted.rds"))

print("We have saved genetic data, we are now removing withdrawal participants")
### Remove withdrawal participants

# Read eid data
withdrawn_eid<-read.csv(paste0(data,"w19266_20200204.csv"))[,1]
ukb_eid<-rownames(genodata)

ukb_eid_clean<-unique(setdiff(ukb_eid,withdrawn_eid))

# Number of participants removed
length(ukb_eid)-length(ukb_eid_clean)

genodata_clean<-genodata[ukb_eid_clean,]
dim(genodata_clean)

### Replace

x<-genodata_clean
x<-sapply(x,as.character)

for (i in seq(1,78)){
  
  x[,i][x[,i] == "03"] <-2
  x[,i][x[,i] == "02"] <-1
  x[,i][x[,i] == "01"] <-0
  x[,i][x[,i] == "00"] <- NA
  
}

# Make df then rename rows
x<-as.data.frame(x)
rownames(x)<-rownames(genodata_clean)

# Check if matching
dim(x)==dim(genodata_clean)

print("Just saving the converted genetic data")
# Save the converted genetic data
saveRDS(x, paste0(res,"genetic_data_converted_all.rds"))
toc()

#----------------- 2c IBD Beta correction -----------------#
# Original code is comparing_alleles.R
# This takes about 2 seconds
print("Starting Step 2c")
tic("Step 2c")

#get relevant columns from both datasets
ukb<-read_bim('import_snps_results_all/ukb_imp_merged.bim')
ukb<-ukb[,c('id','ref','alt')]

#ibd
ibd<-read.table("iv_maf_ibd.txt", header = TRUE, sep = ",", dec = ".")
ibd<-ibd[,c('rsid', 'beta','a1','a2')]
colnames(ibd)<-c('id', 'beta','a1','a2')

#uc
uc<-read.table("iv_maf_uc.txt", header = TRUE, sep = ",", dec = ".")
uc<-uc[,c('rsid', 'beta','a1','a2')]
colnames(uc)<-c('id', 'beta','a1','a2')

#cd
cd<-read.table("iv_maf_cd.txt", header = TRUE, sep = ",", dec = ".")
cd<-cd[,c('rsid', 'beta','a1','a2')]
colnames(cd)<-c('id', 'beta','a1','a2')

print("We've read in all the data - just merging the datasets now :)")
#merge datasets
full_ibd<-merge(ibd, ukb, by='id')
full_uc<-merge(uc, ukb, by='id')
full_cd<-merge(cd, ukb, by='id')

#check whether alleles are the same
full_ibd$check<-(full_ibd$a1==full_ibd$ref & full_ibd$a2== full_ibd$alt)
full_uc$check<-(full_uc$a1==full_uc$ref & full_uc$a2== full_uc$alt)
full_cd$check<-(full_cd$a1==full_cd$ref & full_cd$a2== full_cd$alt)

#if not, check whether they are just switched
full_ibd$check2<-ifelse(full_ibd$check==TRUE, TRUE, ifelse((full_ibd$a2==full_ibd$ref & full_ibd$a1== full_ibd$alt),TRUE,FALSE))
sum(full_ibd$check2) #all are either same or just switched
full_uc$check2<-ifelse(full_uc$check==TRUE, TRUE, ifelse((full_uc$a2==full_uc$ref & full_uc$a1== full_uc$alt),TRUE,FALSE))
sum(full_uc$check2)
full_cd$check2<-ifelse(full_cd$check==TRUE, TRUE, ifelse((full_cd$a2==full_cd$ref & full_cd$a1== full_cd$alt),TRUE,FALSE))
sum(full_cd$check2)

print("Inverting the sateb!")
#invert beta to make it suitable for calculation of PRS with UKB ref alleles
full_ibd$newbeta <-ifelse (full_ibd$check==TRUE, -1*full_ibd$beta, full_ibd$beta)
full_uc$newbeta <-ifelse (full_uc$check==TRUE, -1*full_uc$beta, full_uc$beta)
full_cd$newbeta <-ifelse (full_cd$check==TRUE, -1*full_cd$beta, full_cd$beta)

print("...and just saving the output")
#save output
saveRDS(full_ibd, file='corrected_betas_ibd.rds')
saveRDS(full_uc, file='corrected_betas_uc.rds')
saveRDS(full_cd, file='corrected_betas_cd.rds')

toc()
#----------------- 2d Calculate PRS -----------------#
# Original code is calculate_prs_parallelised_all.R
# This takes about 4 mins to complete

print("Starting step 2d")
tic("Step 2d")

#load beta coefficients for calculation
betas_ibd<-readRDS('corrected_betas_ibd.rds')
betas_ibd$id<-as.character(betas_ibd$id)
rownames(betas_ibd)<-betas_ibd$id

betas_uc<-readRDS('corrected_betas_uc.rds')
betas_uc$id<-as.character(betas_uc$id)
rownames(betas_uc)<-betas_uc$id

betas_cd<-readRDS('corrected_betas_cd.rds')
betas_cd$id<-as.character(betas_cd$id)
rownames(betas_cd)<-betas_cd$id


#load dataset
data<-readRDS('import_snps_results_all/genetic_data_converted_all.rds')

print("Loaded in the data, now saving some col names and making it numeric")

#make small dummy data
#data<-data[1:100,]
#data[2,2]<-NA
#data[6,1:3]<-NA

#save colnames, no of individuals and snps for each disease
col<-colnames(data)
row<-nrow(data)
snps_ibd<-betas_ibd$id
snps_uc<-betas_uc$id
snps_cd<-betas_cd$id


#make data numeric
for (i in col) {
  data[,i]<-as.numeric(data[,i])
}

print("Launching cluster!")
#start cluster
no_cores=min(detectCores(), nchunks)
cl <- makeCluster(no_cores, type="FORK")

prs_ibd<-parApply(cl=cl, X=data, MARGIN = 1,FUN=function(x){
  prs=0
  no_na=0
  #loop over all columns, multiply number of alleles with betavalue for that SNP, just for non-NA SNPS
  for(snp in snps_ibd ){
    if (!is.na(x[snp])){
      prs=prs+ (x[snp]*betas_ibd[snp,'newbeta'])
      no_na=no_na+1
    } 
  }
  #return PRS if <50% snps missing, divide by number of measured SNPs to make PRS comparable between indiv with different number of measured SNPs
  if(no_na>=0.5*length(snps_ibd)){
    return(prs/no_na)
  }  else {
    return(NA)
  }
})

print("Making good progress in the cluster")

prs_uc<-parApply(cl=cl, X=data, MARGIN = 1,FUN=function(x){
  prs=0
  no_na=0
  #loop over all columns, multiply number of alleles with betavalue for that SNP, just for non-NA SNPS
  for(snp in snps_uc ){
    if (!is.na(x[snp])){
      prs=prs+ (x[snp]*betas_uc[snp,'newbeta'])
      no_na=no_na+1
    } 
  }
  #return PRS if <50% snps missing, divide by number of measured SNPs to make PRS comparable between indiv with different number of measured SNPs
  if(no_na>=0.5*length(snps_uc)){
    return(prs/no_na)
  }  else {
    return(NA)
  }
})

print("One more function left to do in the cluster")

prs_cd<-parApply(cl=cl, X=data, MARGIN = 1,FUN=function(x){
  prs=0
  no_na=0
  #loop over all columns, multiply number of alleles with betavalue for that SNP, just for non-NA SNPS
  for(snp in snps_cd ){
    if (!is.na(x[snp])){
      prs=prs+ (x[snp]*betas_cd[snp,'newbeta'])
      no_na=no_na+1
    } 
  }
  #return PRS if <50% snps missing, divide by number of measured SNPs to make PRS comparable between indiv with different number of measured SNPs
  if(no_na>=0.5*length(snps_cd)){
    return(prs/no_na)
  }  else {
    return(NA)
  }
})


#stop cluster
stopCluster(cl)

print("Cluster has landed")

#add PRS to dataframe
data$PRS_ibd<-prs_ibd
data$PRS_uc<-prs_uc
data$PRS_cd<-prs_cd
data$eid<-rownames(data)

#create final part to be exported
prs_final_all<-data[,c('eid','PRS_ibd', 'PRS_uc', 'PRS_cd')]

print("Just need to save it now")

saveRDS(prs_final_all, file='PRS_all.rds')

toc()

print("Completed Step 2 (Genetic Data Creation), the file PRS_all.rds should be in your data folder at full_scripts/data/PRS_all.rds")

toc()

