#----------------- 1: SNP List Creation -----------------#

# This script is complete and takes about 90 mins to run

# This script performs the actions of SNP List Creation, defined by this flowchart: https://whimsical.com/tds-r-scripts-and-data-flow-VmAm6BzY1jUML2a32569t2
# It is designed to run on the HPC server

# It requires
# > IBD Genetics Consortium 3x ".assoc" files
# > UKB Genotyping Data ".bim" and ".fam" files
# > Several packages
library(genio)
library(parallel)
library(data.table)
library(ieugwasr)
library(dplyr)
library(tictoc)

args=commandArgs(trailingOnly=TRUE)
nchunks=as.numeric(args[1])
path=as.character(args[2])
print(path)
data<-paste0(path,"/data/")
setwd(data)
print("We've set the working directory")

print("Welcome to the start of the Job, we are initiating Step 1 now (SNP List Creation)...")
tic("Step 1")

# It outputs a MAFs for IBD-significant SNPs that are also present in UKB
# and a list of these SNPs

#----------------- 1a Get UKB SNPs -----------------#
# Original code is extract_snps_ukb_parallelised.R by Marie 
# This step works and takes between 34 & 80 mins to run
print("Starting 1a")
tic("Step 1a")

#create List of all .bim files
List=NULL
for (chr in 1:22) {
  path<-paste0('/rds/general/project/uk-biobank-2018/live/reference/sdata_12032018/ukb_imp_chr',chr,'.bim')
  List<-rbind(List, path)
}

#make cluster

print("Starting Cluster job")

no_cores=min(detectCores(), nchunks)
cl <- makeCluster(no_cores, type="FORK")
#extract rs_IDs from all .bim files
snps = parLapply(cl=cl, X=List, function(x) {
  bim <- read_bim(x)
  snps_chr<-bim$id
  return(snps_chr)})

stopCluster(cl)

print("Cluster job finished")
print("There's a long wait here, so get some tea and have a nice rest while you wait :)")

#convert list into dataframe
snps_ukb<- data.frame(matrix(unlist(snps), ncol=1, byrow=TRUE))
colnames(snps_ukb)<-'rs_ID'

#save data
#saveRDS(snps_ukb, file='/rds/general/project/hda_students_data/live/Group4/General/full_scripts/data/snps_imp_ukb.rds')
toc()
#----------------- 1b Filter UKB matched SNPs before clumping -----------------#
# Original code is ibdgwas_explore.R by Xin
# This takes about 5 mins to run
print("Starting 1b")
tic("Step 1b")

### Load data ----------------

# Define column index for extraction
colids<-c(2,9,10,11,4,5,6,7)  # "SNP","OR","SE","P","A1","A2","FRQ*"

diseases<-c("ibd","cd","uc")

# <---- read IBD imputed summary level data below
ibd_summary=data.frame(fread("data/EUR.IBD.gwas_info03_filtered.assoc",
                             select=colids))

cd_summary=data.frame(fread("data/EUR.CD.gwas_info03_filtered.assoc",
                            select=colids))

uc_summary=data.frame(fread("data/EUR.UC.gwas_info03_filtered.assoc",
                            select=colids))

# Rename columns for ld_clump()
names(ibd_summary)<-c("rsid","beta","se","pval","a1","a2","freq_a","freq_u")
names(cd_summary)<-c("rsid","beta","se","pval","a1","a2","freq_a","freq_u")
names(uc_summary)<-c("rsid","beta","se","pval","a1","a2","freq_a","freq_u")

# Extract ukb SNPs
#snps_ukb <- readRDS("data/snps_imp_ukb.rds")

### Process summary level data ----------

for (i in 1:length(diseases)){
  
  # Extract disease type
  disease<-diseases[i]
  
  # Parse GWAS summary
  gwas_summary<-eval(parse(text=paste0(disease,"_summary")))
  
  # Rename columns for ld_clump()
  names(gwas_summary)<-c("rsid","beta","se","pval","a1","a2","freq_a","freq_u")
  
  # Filter GWS SNPs for IV
  iv<-gwas_summary[gwas_summary$pval<5*10^-8,]
  
  # Find snps in both ukb and ibd genetics
  snps_mutual<-intersect(snps_ukb[,1],iv$rsid)
  
  # Trim the iv using snps mutually existent in ukb and ibd genetics
  iv_mutual<-iv[iv$rsid %in% snps_mutual,]
  
  # Clump the IV
  iv_clump<-ld_clump(iv_mutual, clump_r2 = 0.01)
  
  # Assign clumped IV to object
  assign(paste0("iv_clump_", disease),iv_clump)
  
  # Filter rare variants
  iv_maf<-filter(iv_clump,freq_a>=0.01&freq_u>=0.01)
  
  iv_maf$beta<-log(iv_maf$beta)
  
  # Assign IV post rare variant removal to object
  assign(paste0("iv_maf_", disease),iv_maf)
  
  # Make export path
  path<-paste0("data/","iv_maf_", disease,".txt")
  
  # Export IV
  write.csv(iv_maf,path)
  
  # Print the IV dimension
  print(paste0("IBD genetics has ", dim(iv)[1]," SNPs significantly associated with ", disease))
  print(paste0("IBD genetics has ", dim(iv_mutual)[1]," SNPs of the above also exist in UKB"))
  print(paste0(disease, ": clumped ", dim(iv_mutual)[1]-dim(iv_clump)[1]," SNPs"))
  print(paste0(disease, ": found ", dim(iv_clump)[1]-dim(iv_maf)[1]," rare SNPs"))
  print(paste0(disease, ", returned ",dim(iv_maf)[1]," SNPs"))
  
}

toc()

#----------------- 1c Create SNP List -----------------#
# Original code is create_list_of_all_snps.R by Marie
# This takes about 1 second to run
print("Starting 1c")
tic("Step 1c")

uc<-read.delim('data/iv_maf_uc.txt', sep = ",", dec = ".")
cd<-read.delim('data/iv_maf_cd.txt', sep = ",", dec = ".")
ibd<-read.delim('data/iv_maf_ibd.txt', sep = ",", dec = ".")

two_rsid<-full_join(uc,cd, by='rsid')
all_rsid<-full_join(two_rsid,ibd, by='rsid')

snps<-all_rsid$rsid
write.table(snps,"data/snps_imp_list_all.txt", row.names = FALSE, col.names=FALSE,quote=FALSE)

toc()

print("Step 1 is Complete, you should have your output files in your data folder")
print("Now initiating Step 2, Genetic Data Creation")
toc()
