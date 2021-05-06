#----------------- 5PLUS: Analysis & Visualisation -----------------#

# This script is complete and takes about XXXX mins to run

# This script performs the actions of Analysis & Visualisation, defined by this flowchart: https://whimsical.com/tds-r-scripts-and-data-flow-VmAm6BzY1jUML2a32569t2
# It is designed to run on the HPC server

print("Now launching RScript 5PLUS!")

# It requires
# > results from 5_Analysis_&_Visualisation.R
# > Several packages
library(tidyverse)#installed
library(tictoc)#installed
library(regclass) #installed
library(car) #installed
library(imputeMissings) #for imputation #installed
library(oem) # installed
library(parallel)#installed
library(matchmaker) #installed


args=commandArgs(trailingOnly=TRUE)
nchunks=as.numeric(args[1])
path=as.character(args[2])
print(path)
datapath<-paste0(path,"/data/")
analysis<-paste0(path,"/analysis/")
setwd(datapath)
print("We've set the working directory")

print("Welcome to the start of the Job, we are initiating the end of Step 5 now!!")
tic("Step 5plus")

# It outputs analysis and visualisation in the folder "analysis"


#----------------- 5g LASSO -----------------#
# Original code is lasso_oem.R by Marie 
print("Starting step 5g")
tic("Step 5g")

#read data already with dummy variables
newdata<-readRDS('matrixLasso.rds')
#create X and y data
drop<-c('colon1','ibd1', '(Intercept)')
X<-newdata[, !colnames(newdata) %in% drop]
y_c<-newdata[,'colon1'] #for gLASSO on colon cancer
y_ibd<-newdata[,'ibd1'] #for gLASSO on ibd

#X<-X[1:100,]
#y_c<-y_c[1:100]
#y_c[5]<-1
#y_c[34]<-1
#y_c[78]<-1
#y_c[79]<-1
#y_c[80]<-1
#y_ibd<-y_ibd[1:100]

print("Data loaded, now creating corresponding group file")
#create corresponding group file
cols<-colnames(X) 
groups<-1 
g=1
for (i in 2:(length(cols))) { #loop over all colnames starting from position2
  if (gsub("[0-9]", "", x=cols[i]) == gsub("[0-9]", "", x=cols[i-1])){ #check whether name of column is same as previous colname (without indices)
  }else{
    g=g+1 #if not same name: new group
  }
  groups<-c(groups,g) #add group number to group list
}
print("Now fitting group Lasso incl. Stability analysis")
#fit group lasso incl. stability analysis
gLassoSub = function(k , Xdata, Ydata, Groups, Nfolds=5) {
  set.seed(k)
  s = sample(nrow(Xdata), size = 0.8 * nrow(Xdata)) #take random subsample of 80%
  Xsub = Xdata[s, ]
  Ysub = Ydata[s]
  #cross-validation
  cvfit1 <- cv.oem(x = Xsub, y = Ysub, penalty = c("grp.lasso"), type.measure = "deviance", groups = Groups, nfolds = Nfolds)
  #fitting model with lambda1se
  fit <- oem(x = Xsub, y = Ysub, penalty = c("grp.lasso"), family = "binomial", groups = Groups, lambda = cvfit1$lambda.1se.models)
  return(fit$beta$grp.lasso) 
}
print("Launching parallelisation now")
t0=Sys.time()
#start cluster for parallelization
no_cores=min(detectCores(), nchunks)
cl <- makeCluster(no_cores, type="FORK")

#run stability analysis
niter = 1000 # iterations for stability analysis
lasso.stab_ibd = parSapply(cl=cl, 1:niter, FUN = gLassoSub, Xdata = X,   Ydata = y_ibd, Groups=groups, Nfolds=5)
print("Half way through the parallelisation")
lasso.stab_c = parSapply(cl=cl, 1:niter, FUN = gLassoSub, Xdata = X,   Ydata = y_c, Groups=groups, Nfolds=5)

#stop cluster
stopCluster(cl)

t1=Sys.time()
print("End of parallelisation")
computation_time=difftime(t1,t0, units="secs")
print(computation_time)

rownames(lasso.stab_ibd)<-c('Intercept',cols) #add variable names to dataset
rownames(lasso.stab_c)<-c('Intercept',cols) #add variable names to dataset

print("Now exploring results")
#------------------explore results-------------------------#
countNonzero<-function(x,dat){
  sum(dat[x,]!=0)
}

#summ up selections across all runs
sums_colon<-data.frame(sapply(1:237, FUN=countNonzero, dat=lasso.stab_c))
sums_ibd<-data.frame(sapply(1:237, FUN=countNonzero, dat=lasso.stab_ibd))
rownames(sums_colon)<- rownames(lasso.stab_c)
rownames(sums_ibd)<- rownames(lasso.stab_ibd)


result_colon<-data.frame(count=sums_colon[1,])
rows<-rownames(sums_colon)
s=2
for (i in 2:(length(rownames(sums_colon)))) { 
  if (gsub("[0-9]", "", x=rows[i]) == gsub("[0-9]", "", x=rows[i-1])){  #check whether name of row is same as previous rowname (without indices)
    check<-sums_colon[i,]==sums_colon[i-1,] #check whether selection count was the same for those
    #print(check)
  }else{
    result_colon[s,]<-sums_colon[i,]
    rownames(result_colon)[s]<-gsub("[0-9]", "", x=rownames(sums_colon)[i]) #delete numbers at the end of names
    s=s+1 #next row
  }
}

result_ibd<-data.frame(count=sums_ibd[1,])
rows<-rownames(sums_ibd)
s=2
for (i in 2:(length(rownames(sums_ibd)))) { 
  if (gsub("[0-9]", "", x=rows[i]) == gsub("[0-9]", "", x=rows[i-1])){  #check whether name of row is same as previous rowname (without indices)
    check<-sums_ibd[i,]==sums_ibd[i-1,] #check whether selection count was the same for those
    #print(check)
  }else{
    result_ibd[s,]<-sums_ibd[i,]
    rownames(result_ibd)[s]<- gsub("[0-9]", "", x=rownames(sums_ibd)[i]) #delete numbers at the end of names
    s=s+1 #next row
  }
}

#filter for nonzero variables
result_ibd_nzero<-filter(result_ibd, count != 0)
ibd_vars<-rownames(result_ibd_nzero)[-1]
result_colon_nzero<-filter(result_colon, count != 0)
colon_vars<-rownames(result_colon_nzero)[-1]

#extract intersect
intersect<- ibd_vars[ibd_vars%in% colon_vars]
print("Just saving RDS now")
#save data
saveRDS(intersect, file=paste0(datapath,'Lasso_intersect.rds'))

toc()

#----------------- 5h Attenuation Analysis -----------------#
# Original code is logreg_colon.R by Marie
print("Starting step 5h")
tic("Step 5h")

#set path

#read data and specify factors/numerical variables
raw_data<-readRDS('ukb_hes_everything.rds')
continuous<-c("eid", "PRS_ibd", "PRS_uc", "PRS_cd", "age", "bmi", "duration_walks", "days_moderate_activity_per_week", "days_vigorous_activity_per_week",
              "days_walk_per_week", "sleep_duration", "exposure_tobacco_smoke_home", "exposure_tobacco_smoke_outside", "cooked_veg", "bread_intake", "salad_raw_veg",
              "cereal", "tea", "coffee", "water", "redwine_weekly", "champagnewhitewine_weekly", "beercide_weekly", "spirits_weekly", "fortifiedwine_weekly",
              "age_secual_intercourse", "lifetime_sexual_partners", "albumin", "alanine_aminotransferase", "alkaline_phosphatase", "apolipoprotein_A", "apolipoprotein_B", "aspartate_aminotransferase",
              "C-reactive_protein", "calcium", "cholesterol", "creatinine", "cystatin_C", "direct_bilirubin", "gamma_glutamyltransferase", "glucose", "glycated_haemoglobin(HbA1c)",
              "hdl_cholesterol", "igf-1", "ldl_direct", "lipoprotein_A", "phosphate", "shbg", "testosterone", "total_bilirubin", "total_protein", "triglycerides", "urate",
              "urea", "vitamin_D", "duration_moderate_activity", "fresh_fruit", "dried_fruit")
factors<-colnames(raw_data) [! colnames(raw_data) %in% continuous] #everything else
raw_data[,factors] <- lapply(raw_data[,factors] , factor)
raw_data[,continuous] <- lapply(raw_data[,continuous] , as.numeric)

#impute missing values with median/mode
imp_data<-imputeMissings::impute(raw_data)

###unadjusted model
fit_unadjusted<-glm(colon~ibd, family = binomial(link=logit), data=imp_data)
#summary(fit_unadjusted)
#exp(coef(summary(fit_unadjusted)))

#in case variable names were renamed during analyses:, correcting them so they fit the imputed dataset version
corrections <- data.frame(
  bad = c(".regex handsfree_device_last([0-9]*|.)mo",".regex alcohol_vs_([0-9]*|.)yrs",".regex (X*|5|.)asa",".regex C(.*|-)reactive_protein",".regex glycated_haemoglobin(\\(|.*)HbA(1*|.)c(\\(|.*)",".regex igf(-*|.|[0-9])"),
  good = c("handsfree_device_last3mo", "alcohol_vs_10yrs", "X5asa", "C.reactive_protein", "glycated_haemoglobin.HbA1c.", "igf.1") ,
  stringsAsFactors = FALSE
)
print("Reading in LASSO selected variables")
###only lasso selected variables
vars_lasso <- readRDS("Lasso_intersect.rds")
vars_lasso<-c('colon','ibd', match_vec(vars_lasso, corrections)) #correcting changed var names and adding colon and ibd
#vars_lasso<-c('colon','ibd','major_diet_changes','loperamid','alkaline_phosphatase','cystatin_C','unspecific_benign_neoplasm_of_colon','FH_MNDO')
data_lasso<-imp_data[,(names(imp_data) %in% vars_lasso)]
fit_adjusted_lasso<-glm(colon~., family = binomial(link=logit), data=data_lasso)
summary(fit_adjusted_lasso)
exp(coef(summary(fit_adjusted_lasso)))

print("Reading in random forest selected variables")
###random forest selected variables
rf_pred <- readRDS("ukb_ML_randomforest_jointpredictors_perm.rds")
rf_pred[] <- lapply(rf_pred, as.character)
vars_rf<-unlist(c('colon','ibd',rf_pred$ibd_gini_predictors))
vars_rf<-vars_rf[!vars_rf %in% c('ldl_direct','cholesterol')] #issue with collinearity!!
vars_rf<-match_vec( vars_rf, corrections) #correcting changed var names
data_rf<-imp_data[,(names(imp_data) %in% vars_rf)]
fit_adjusted_rf<-glm(colon~., family = binomial(link=logit), data=data_rf)
#summary(fit_adjusted_rf)
#exp(coef(summary(fit_adjusted_rf)))
#vif(fit_adjusted_rf)

print("Reading in sPLS selected variables")
###sPLS selected variables
#read in Nas' file? somehow not possible
vars_pls<-c('colon','ibd',"age","smoking_status","household_income","days_vigorous_activity_per_week","usual_walking_pace","nap_during_day", 
            "sleeplessness","ever_smoked","past_tobacco_smoking","bread_type","cereal_type","coffee_type","major_diet_changes" ,"loperamid",
            "codeine_phosphate","albumin","alkaline_phosphatase","C.reactive_protein","calcium","cholesterol","creatinine","cystatin_C",
            "gamma_glutamyltransferase","hdl_cholesterol","ldl_direct","unspecific_benign_neoplasm_of_colon","FH_MNDO","Depressive_episode",
            "depression","diabtes_merge", "never_eat_egg_dairy_wheat_sugar", "alcohol_intake_freq") 
data_pls<-imp_data[,(names(imp_data) %in% vars_pls)]
fit_adjusted_pls<-glm(colon~., family = binomial(link=logit), data=data_pls)
#summary(fit_adjusted_pls)
#exp(coef(summary(fit_adjusted_pls)))


print("Writing to file")
# Write to file
sink(file=paste0(analysis,"Attenuation_logistic_regressions_results.txt"), type="output", append=TRUE)
print('Unadjusted model:')
summary(fit_unadjusted)
cat(paste(" ", " ", " "," ", " ", "", sep="\n"))
print('Adjusted model, sPLS selected variables:')
summary(fit_adjusted_pls)
cat(paste(" ", " ", " "," ", " ", "", sep="\n"))
print('Adjusted model, Lasso selected variables:')
summary(fit_adjusted_lasso)
cat(paste(" ", " ", " "," ", " ", "", sep="\n"))
print('Adjusted model, Random Forest selected variables:')
summary(fit_adjusted_rf)
sink() # stop writing to file


toc()

print("Step 5 (fully) complete!")
toc()