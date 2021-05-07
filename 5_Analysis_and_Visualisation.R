#----------------- 5: Analysis & Visualisation -----------------#

# This script is complete and takes about 28+RF mins to run

# This script performs the actions of Analysis & Visualisation
# It is designed to run on the HPC server

# It requires
# > ukb_hes_everything.rds
# > Several packages
library(tidyverse)
library(tictoc)
library(gtsummary)
library(gt)
library(data.table)
library(patchwork)
library(ggplot2)
library(plotROC)
library(rcompanion) #for nagelkerke 
library(imputeMissings)
library(sgPLS)
library(VennDiagram)
library(pheatmap)
library(utils)
library(ROSE) 
library(gsubfn) 
library(randomForest)
library(gridExtra) 
library(Hmisc)
library(RColorBrewer)

args=commandArgs(trailingOnly=TRUE)
nchunks=as.numeric(args[1])
path=as.character(args[2])
print(path)
datapath<-paste0(path,"/data/")
analysis<-paste0(path,"/analysis/")
setwd(datapath)
print("We've set the working directory")

print("Welcome to the start of the Job, we are initiating Step 5 now (Analysis & Visualisation)...")
tic("Step 5")

# It outputs analysis and visualisation in the folder "analysis"

#----------------- 5a Table 1 -----------------#
# Original code is table1.R by Xin
# This step works and takes between less than 1 min to run
print("Starting 5a")
tic("Step 5a")

# Load the data
ukb_hes_everything<-readRDS('ukb_hes_everything.rds')
raw_data <- ukb_hes_everything
df <- ukb_hes_everything

print("tidying data")
# Tidy data
df$gender<-as.factor(ifelse(df$gender==1,"male","female"))
df$smoking_status<-ifelse(is.na(df$smoking_status),"prefer not to say",
                            ifelse(df$smoking_status==2,"current",
                                   ifelse(df$smoking_status==1,"previous",
                                          "never")))
df$smoking_status<-factor(df$smoking_status, 
                            levels=c("current","previous","never","prefer not to say"))
df$colon<-as.factor(ifelse(df$colon==1,"Case","Control"))
print("standardising PRS scores")
# Standardise PRS scores
df[,(names(df) %in% c('PRS_ibd','PRS_cd','PRS_uc'))] <- scale( df[,(names(df) %in% c('PRS_ibd','PRS_cd','PRS_uc'))] )

# Specify columns to include
cname_inc<-c("colon","ibd","uc","cd","age","gender","bmi","smoking_status","PRS_ibd","PRS_uc","PRS_cd")
# Produce table1 <------------ If troubles installing packages to run the following code, use "tableone"
table1<-df %>%
  dplyr::select(cname_inc) %>%
  tbl_summary(
    by = colon,
    statistic = list(all_continuous() ~ "{mean} ({sd})",
                     all_categorical() ~ "{n} ({p}%)"),
    missing = "no",
    label = list(
      smoking_status ~ "smoking status",
      PRS_ibd ~ "prs ibd", PRS_cd ~ "prs cd", PRS_uc ~ "prs uc"
    ),
    digits = all_continuous() ~ 2,
    missing_text = "(Missing)"
  ) %>%
  add_overall() %>%
  modify_header(label = "**Variable**") %>%
  bold_labels() 

# Export table as image
gtsave(as_gt(table1), path = analysis, "table1.png")


toc()

#----------------- 5b Logistic Regression -----------------#
# Original code is logistic_regression.R 
# This step works and takes less than 1 min to run
print("Starting 5b")
tic("Step 5b")

#read data and specify factors/numerical variables
# raw_data already read in previously
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

#use only participants with genetic data available
data<-raw_data[!is.na(raw_data$PRS_ibd),]

#standardize PRS
data[,(names(data) %in% c('PRS_ibd','PRS_cd','PRS_uc'))] <- scale( data[,(names(data) %in% c('PRS_ibd','PRS_cd','PRS_uc'))] )

#-------------------------------logistic regression colon cancer~PRS_ibd---------------#
#check asumption: linear relationship of logodds(colon cancer) and PRS
data$PRS_ibd_group <- as.factor(cut2(data$PRS_ibd, g=50))
dict <-aggregate(data$PRS_ibd, by= list(data$PRS_ibd_group), FUN=mean)
map = setNames(dict$x, dict$Group.1)
data$PRS_ibd_mean <- map[unlist(data$PRS_ibd_group)]
cc_by_prs<-as.data.frame.matrix(table(data$PRS_ibd_mean, data$colon))
cc_by_prs$logodds<-log(cc_by_prs[,'1']/cc_by_prs[,'0'])
cc_by_prs$PRS<-as.numeric(rownames(cc_by_prs))
cc_scatter<-ggplot(data=cc_by_prs,aes(x=PRS,y=logodds))+
  geom_point()+
  theme_classic()+
  xlab('IBD PRS')+
  ylab('log(odds) of having Colon Cancer')
#plot(cc_scatter)

#density plot cc cases vs. non cases
cc_dens <- ggplot(data, aes(x=PRS_ibd, colour= colon, fill=colon))+
  geom_density(aes(y=..density..),alpha=0.2, size = 1) +
  scale_fill_manual(values=c("#69b3a2", "#404080"),labels = c("Controls", "Colon Cancer Cases"), name="") +
  scale_colour_manual(values = c("#69b3a2", "#404080"),labels = c("Controls", "Colon Cancer Cases"), name="") +
  theme_classic()+
  theme(legend.position = c(.0, 1),
        legend.justification = c("left", "top"),
        legend.box.just = "left")+
  xlab('IBD PRS')+
  ylab('density')
#plot(cc_dens)

patchwork<-(cc_scatter + cc_dens)
patchwork_annot<-patchwork + plot_annotation(tag_levels = list(c('a)','b)')))

# save as jpeg
jpeg(paste0(analysis,"MR_plot.jpg"), width = 700, height = 350)
patchwork_annot
dev.off()

#fit logistic regression model with PRS_ibd to predict colon cancer
fit<-glm(colon~PRS_ibd, family = binomial(link=logit), data=data)
summary(fit) 
exp(coef(summary(fit)))

#------------------look at other cancer outcomes---------------------------------#
fit_small_intestine<-glm(small_intestine~PRS_ibd, family = binomial(link=logit), data=data)
summary(fit_small_intestine)   
exp(coef(summary(fit_small_intestine)))

fit_rectosigmoid<-glm(rectosigmoid~PRS_ibd, family = binomial(link=logit), data=data)
summary(fit_rectosigmoid)   
exp(coef(summary(fit_rectosigmoid)))

fit_rectal<-glm(rectal~PRS_ibd, family = binomial(link=logit), data=data)
summary(fit_rectal) 
exp(coef(summary(fit_rectal)))

fit_anal<-glm(anal~PRS_ibd, family = binomial(link=logit), data=data)
summary(fit_anal) 
exp(coef(summary(fit_anal)))

fit_liver<-glm(liver~PRS_ibd, family = binomial(link=logit), data=data)
summary(fit_liver)
exp(coef(summary(fit_liver)))


#---------------look specifically at small_intestine~PRS_cd and rectosigmoid~PRD_uc-----------#

#small intestine cancer~PRS_cd
fit_small_intestine_cd<-glm(small_intestine~PRS_cd, family = binomial(link=logit), data=data)
summary(fit_small_intestine_cd)
exp(coef(summary(fit_small_intestine_cd)))

#rectosigmoid cancer~PRS_uc
fit_rectosigmoid_uc<-glm(rectosigmoid~PRS_uc, family = binomial(link=logit), data=data)
summary(fit_rectosigmoid_uc) 

#----#

##small intestine cancer~PRS_uc
fit_small_intestine_uc<-glm(small_intestine~PRS_uc, family = binomial(link=logit), data=data)
summary(fit_small_intestine_uc)
exp(coef(summary(fit_small_intestine_uc)))

#rectosigmoid cancer~PRS_cd
fit_rectosigmoid_cd<-glm(rectosigmoid~PRS_cd, family = binomial(link=logit), data=data)
summary(fit_rectosigmoid_cd)


# Write logistic regression results to file
sink(file=paste0(analysis,"MR_logistic_regressions_results.txt"), type="output", append=TRUE)
print('MR analyses: PRS-IBD related to different digestive cancers')
cat(paste(" ", "", sep="\n"))
print('colon cancer:')
summary(fit)
cat(paste(" ", " ", " "," ", " ", "", sep="\n"))
print('small intestine cancer:')
summary(fit_small_intestine) 
cat(paste(" ", " ", " "," ", " ", "", sep="\n"))
print('rectosigmoid cancer:')
summary(fit_rectosigmoid)
cat(paste(" ", " ", " "," ", " ", "", sep="\n"))
print('rectal cancer:')
summary(fit_rectal) 
cat(paste(" ", " ", " "," ", " ", "", sep="\n"))
print('anal cancer:')
summary(fit_anal) 
cat(paste(" ", " ", " "," ", " ", "", sep="\n"))
print('liver cancer:')
summary(fit_liver) 
cat(paste(" ", " ", " "," ", " ", "", sep="\n"))
print('MR analyses: PRS-UC and PRS-CD related to different small intestine and rectosigmoid cancer')
cat(paste(" ", "", sep="\n"))
print('small intestine cancer by PRS-CD:')
summary(fit_small_intestine_cd)
cat(paste(" ", " ", " "," ", " ", "", sep="\n"))
print('small intestine cancer by PRS-UC:')
summary(fit_small_intestine_uc)
cat(paste(" ", " ", " "," ", " ", "", sep="\n"))
print('rectosigmoid cancer by PRS-CD:')
summary(fit_rectosigmoid_cd)
cat(paste(" ", " ", " "," ", " ", "", sep="\n"))
print('rectosigmoid cancer by PRS-UC:')
summary(fit_rectosigmoid_uc) 
sink() # stop writing to file
toc()

#----------------- 5c QC of PRS -----------------#
# Original code is visualisation_onlyPRS.R
# This step works and takes less than 1 min to run
print("Starting 5c")
tic("Step 5c")

# data df has already been created

################################## Quality Control of PRS ################################


#-------------------get mean prs by disease status-----------------------#
aggregate(data$PRS_ibd, list(data$ibd), mean)
aggregate(data$PRS_cd, list(data$cd), mean)
aggregate(data$PRS_uc, list(data$uc), mean)
aggregate(data$PRS_ibd, list(data$colon), mean)

mean(data$PRS_ibd)

#############PRS_IBD####################

#logodds of ibd by prs
#make subgroups for plotting 
data$PRS_ibd_group <- as.factor(cut2(data$PRS_ibd, g=50))
dict <-aggregate(data$PRS_ibd, by= list(data$PRS_ibd_group), FUN=mean)
map = setNames(dict$x, dict$Group.1)
data$PRS_ibd_mean <- map[unlist(data$PRS_ibd_group)]

ibd_by_prs<-as.data.frame.matrix(table(data$PRS_ibd_mean, data$ibd))
ibd_by_prs$logodds<-log(ibd_by_prs[,'1']/ibd_by_prs[,'0'])
ibd_by_prs$PRS<-as.numeric(rownames(ibd_by_prs))
ibd1<-ggplot(data=ibd_by_prs,aes(x=PRS,y=logodds))+
  geom_point()+
  theme_classic()+
  xlab('IBD PRS')+
  ylab('log(odds) of having IBD')
#plot(ibd1)

# plot distribution of PRS grouped by ibd status
ibd2 <- ggplot(data, aes(x=PRS_ibd, colour= ibd, fill=ibd))+
  geom_density(aes(y=..density..),alpha=0.2, size = 1) +
  scale_fill_manual(values=c("#69b3a2", "#404080"),labels = c("Controls", "IBD Cases"), name="") +
  scale_colour_manual(values = c("#69b3a2", "#404080"),labels = c("Controls", "IBD Cases"), name="") +
  theme_classic()+
  theme(legend.position = c(.0, 1),
        legend.justification = c("left", "top"),
        legend.box.just = "left")+
  xlab('IBD PRS')+
  ylab('density')
#plot(ibd2)

#plot ROC and calculate AUC
rocdata <- data.frame(D =data$ibd,M = data$PRS_ibd)
rocdata$D<-as.integer(rocdata$D)-1
roc_ibd<-ggplot(rocdata, aes(m = M, d = D)) + 
  geom_roc(n.cuts = 0) + style_roc() + geom_abline(intercept = 0, slope = 1,linetype="dashed", size=.5)

ibd3 <- roc_ibd + annotate("text", x = .8, y = .1, label = paste("AUC = ",round(calc_auc(roc_ibd)[3],3)))
#plot(ibd3)


#############PRS_UC####################
#logodds of uc by PRS
data$PRS_uc_group <- as.numeric(cut2(data$PRS_uc, g=50))
dict <-aggregate(data$PRS_uc, by= list(data$PRS_uc_group), FUN=mean)
map = setNames(dict$x, dict$Group.1)
data$PRS_uc_mean <- map[unlist(data$PRS_uc_group)]

uc_by_prs<-as.data.frame.matrix(table(data$PRS_uc_mean, data$uc))
uc_by_prs$PRS<-as.numeric(rownames(uc_by_prs))
uc_by_prs$logodds<-log(uc_by_prs[,'1']/uc_by_prs[,'0'])
uc1<-ggplot(data=uc_by_prs,aes(x=PRS,y=logodds))+
  geom_point()+
  theme_classic()+
  xlab('UC PRS')+
  ylab('log(odds) of having UC')
#plot(uc1)

# plot distribution of PRS grouped by uc status
uc2 <- data %>%
  ggplot( aes(x=PRS_uc, colour= uc, fill=uc))+
  geom_density(aes(y=..density..),alpha=0.2, size = 1) +
  scale_fill_manual(values=c("#69b3a2", "#404080"),labels = c("Controls", "UC Cases"), name="") +
  scale_colour_manual(values = c("#69b3a2", "#404080"),labels = c("Controls", "UC Cases"), name="") +
  theme_classic()+
  theme(legend.position = c(.0, 1),
        legend.justification = c("left", "top"),
        legend.box.just = "left")+
  xlab('UC PRS')+
  ylab('density')
#plot(uc2)

#plot ROC and calculate AUC
rocdata_uc <- data.frame(D =data$uc,M = data$PRS_uc)
rocdata_uc$D<-as.integer(rocdata_uc$D)-1
roc_uc<-ggplot(rocdata_uc, aes(m = M, d = D)) + 
  geom_roc(n.cuts = 0) + style_roc() + geom_abline(intercept = 0, slope = 1,linetype="dashed", size=.5)

uc3 <- roc_uc + annotate("text", x = .8, y = .1, label = paste("AUC = ",round(calc_auc(roc_uc)[3],3)))
#plot(uc3)

#############PRS_CD####################
#logodds of cd by PRS
data$PRS_cd_group <- as.numeric(cut2(data$PRS_cd, g=50))
dict <-aggregate(data$PRS_cd, by= list(data$PRS_cd_group), FUN=mean)
map = setNames(dict$x, dict$Group.1)
data$PRS_cd_mean <- map[unlist(data$PRS_cd_group)]

cd_by_prs<-as.data.frame.matrix(table(data$PRS_cd_mean, data$cd))
cd_by_prs$PRS<-as.numeric(rownames(cd_by_prs))
cd_by_prs$logodds<-log(cd_by_prs[,'1']/cd_by_prs[,'0'])
cd1<-ggplot(data=cd_by_prs,aes(x=PRS,y=logodds))+
  geom_point()+
  theme_classic()+
  xlab('CD PRS')+
  ylab('log(odds) of having CD')
#plot(cd1)

# plot distribution of PRS grouped by cd status
cd2 <- data %>%
  ggplot( aes(x=PRS_cd, colour= cd, fill=cd))+
  geom_density(aes(y=..density..),alpha=0.2, size = 1) +
  scale_fill_manual(values=c("#69b3a2", "#404080"),labels = c("Controls", "CD Cases"), name="") +
  scale_colour_manual(values = c("#69b3a2", "#404080"),labels = c("Controls", "CD Cases"), name="") +
  theme_classic()+
  theme(legend.position = c(.0, 1),
        legend.justification = c("left", "top"),
        legend.box.just = "left")+
  xlab('CD PRS')+
  ylab('density')
#plot(cd2)

#plot ROC and calculate AUC
rocdata_cd <- data.frame(D =data$cd,M = data$PRS_cd)
rocdata_cd$D<-as.integer(rocdata_cd$D)-1
roc_cd<-ggplot(rocdata_cd, aes(m = M, d = D)) + 
  geom_roc(n.cuts = 0) + style_roc() + geom_abline(intercept = 0, slope = 1,linetype="dashed", size=.5)
cd3 <- roc_cd + annotate("text", x = .8, y = .1, label = paste("AUC = ",round(calc_auc(roc_cd)[3],3)))
#plot(cd3)


###----- plot all PRS QC plots together -----##

patchwork<-(ibd1+ibd2+ibd3)/(uc1 +uc2 +uc3)/(cd1+ cd2+ cd3)
patch_annot<-patchwork + plot_annotation(
  title = 'Evaluation of the Polygenic Risk Scores',
  tag_levels = list(c('a) 1', '2 ','3','b)', ' ','','c)', ' ', '')))

# save as jpeg
jpeg(paste0(analysis,"PRS_validation_plot.jpg"), width = 750, height = 550)
patch_annot
dev.off()

#--------------calculate nagelkerke pseudo R-squared----------------------#

fit_nk_ibd<-glm(ibd~PRS_ibd, family = binomial(link=logit), data=data)
summary(fit_nk_ibd)
nagelkerke(fit_nk_ibd)$Pseudo.R.squared.for.model.vs.null
exp(coef(summary(fit_nk_ibd)))

fit_nk_cd<-glm(cd~PRS_cd, family = binomial(link=logit), data=data)
summary(fit_nk_cd)
nagelkerke(fit_nk_cd)$Pseudo.R.squared.for.model.vs.null

fit_nk_uc<-glm(uc~PRS_uc, family = binomial(link=logit), data=data)
summary(fit_nk_uc)
nagelkerke(fit_nk_uc)$Pseudo.R.squared.for.model.vs.null

# Write to file
sink(file=paste0(analysis,"PRSvalidation_logistic_regressions_results.txt"), type="output", append=TRUE)
print('Logistic regressions between PRS and their related disease and their nagelkerke pseudo R-squared:')
print('IBD by PRS IBD:')
summary(fit_nk_ibd)
nagelkerke(fit_nk_ibd)$Pseudo.R.squared.for.model.vs.null
cat(paste(" ", " ", " "," ", " ", "", sep="\n"))
print('CD by PRS CD:')
summary(fit_nk_cd)
nagelkerke(fit_nk_cd)$Pseudo.R.squared.for.model.vs.null
cat(paste(" ", " ", " "," ", " ", "", sep="\n"))
print('UC by PRS UC:')
summary(fit_nk_uc)
nagelkerke(fit_nk_uc)$Pseudo.R.squared.for.model.vs.null
sink() # stop writing to file

toc()

#----------------- 5d Standardisation and Dummy Variable Creation -----------------#
# Original code is createMatrixLasso.R
# This step works and takes 60s to run
print("Starting 5d")
tic("Step 5d")

#read data and specify factors/numerical variables
raw_data<-ukb_hes_everything
print("ok now we're taking continuous vars")
continuous<-c("eid", "PRS_ibd", "PRS_uc", "PRS_cd", "age", "bmi", "duration_walks", "days_moderate_activity_per_week", "days_vigorous_activity_per_week",
              "days_walk_per_week", "sleep_duration", "exposure_tobacco_smoke_home", "exposure_tobacco_smoke_outside", "cooked_veg", "bread_intake", "salad_raw_veg",
              "cereal", "tea", "coffee", "water", "redwine_weekly", "champagnewhitewine_weekly", "beercide_weekly", "spirits_weekly", "fortifiedwine_weekly",
              "age_secual_intercourse", "lifetime_sexual_partners", "albumin", "alanine_aminotransferase", "alkaline_phosphatase", "apolipoprotein_A", "apolipoprotein_B", "aspartate_aminotransferase",
              "C-reactive_protein", "calcium", "cholesterol", "creatinine", "cystatin_C", "direct_bilirubin", "gamma_glutamyltransferase", "glucose", "glycated_haemoglobin(HbA1c)",
              "hdl_cholesterol", "igf-1", "ldl_direct", "lipoprotein_A", "phosphate", "shbg", "testosterone", "total_bilirubin", "total_protein", "triglycerides", "urate",
              "urea", "vitamin_D", "duration_moderate_activity", "fresh_fruit", "dried_fruit")
factors<-colnames(raw_data) [! colnames(raw_data) %in% continuous] #everything else
print("Applying correct categorisation to variables")
raw_data[,factors] <- lapply(raw_data[,factors] , factor)
raw_data[,continuous] <- lapply(raw_data[,continuous] , as.numeric)
print("Now we're imputing values")
tic("Now imputing values")
#impute missing values with median/mode
imp_data<-imputeMissings::impute(raw_data)
saveRDS(imp_data,"imputed_everything.rds")
toc()

#drop irrelevant cols for variable selection
print("Dropping irrelevant columns")
#imp_data<-readRDS("imputed_everything.rds")
#imp_data<-as.data.frame(imp_data)
#print("re-imported rds")
irrelevant<-c("eid", "PRS_ibd", "PRS_uc", "PRS_cd","cd","uc", "small_intestine", "rectosigmoid", "rectal", "anal", "liver")
relevant_data<-imp_data[,!(names(imp_data) %in% irrelevant)]

#standardize continuous values
print("standardising continuous values")
relevant_data[,(names(relevant_data) %in% continuous)] <- scale(relevant_data[,(names(relevant_data) %in% continuous)] )

#make factors dummy variables
print("making factors dummy vars")
newdata<-model.matrix(~.,relevant_data)

#save data
print("saving data")
saveRDS(newdata, file='matrixLasso.rds')

toc()

#----------------- 5e sPLS -----------------#
# Original code is sPLS.R 
# This step works and takes 24 mins to run
print("Starting 5e")
tic("Step 5e")

source("pls_functions.R")

# Define Path 
df<- imp_data #readRDS("imputed_everything.rds")

# Using downsampling to fix calibration issue
downsample_result <- ovun.sample(colon~., data=df, method="under", N=7000, seed=2018)
downsample_result2 <- ovun.sample(ibd~., data=df, method="under", N=14000, seed=2018)
downsample_df<-downsample_result$data
downsample_df2<-downsample_result2$data
#dim(downsample_df)


downsample_m <- data.matrix(downsample_df)
downsample_m2 <- data.matrix(downsample_df2)

#prepare downsampled data for PLS: X & y
drop<-c('colon','ibd', 'eid','PRS_ibd', 'PRS_uc', 'PRS_cd','cd','uc', "small_intestine", "rectosigmoid", "rectal", "anal", "liver")
X_d<-downsample_m[, !colnames(downsample_m) %in% drop]
X_d2 <-downsample_m2[, !colnames(downsample_m2) %in% drop]
y_dc<-downsample_m[,'colon'] 
y_dibd<-downsample_m2[,'ibd'] 
y_dcfactor <-as.factor(y_dc)
y_dibdfactor <- as.factor(y_dibd)

table(downsample_df$colon)
table(downsample_df2$ibd)



# perform 5-fold cross-validation to calibrate the number of variables selected by sPLS- DA.
set.seed(1)
res_splsda_c=CalibratesPLSDA(dataX=X_d, dataY=y_dcfactor, ncomp=1, Nrepeat=10)
# Plot the calibration resulr for colon cancer
#pdf(file = "../Output/calibspls_c.pdf", width = 9,height = 6)
PlotCalib(res=res_splsda_c)
dev.off()


set.seed(1)
res_splsda_ibd=CalibratesPLSDA(dataX=X_d2, dataY=y_dibdfactor, ncomp=1, Nrepeat=10)
# Plot the calibration result for IBD
#pdf(file = "../Output/calibspls_ibd.pdf", width = 9,height = 6)
PlotCalib(res=res_splsda_ibd)
dev.off()

# Opitimal number of variables to be selected for Colon cancer: 
# Using sPLS on full data
newdata<-readRDS("matrixLasso.rds")

#prepare data for sPLS: X & y
drop<-c('colon1','ibd1', 'eid','PRS_ibd', 'PRS_uc', 'PRS_cd','cd','uc', "small_intestine", "rectosigmoid", "rectal", "anal", "liver")
X<-newdata[, !colnames(newdata) %in% drop]
y_c<-newdata[,'colon1'] 
y_ibd<-newdata[,'ibd1'] 
y_cfactor <-as.factor(y_c)
y_ibdfactor <- as.factor(y_ibd)
nlevels(y_cfactor)

#Running a calibrated sparse PLS-DA model where keepX is the optimal number of variables selected
MysPLSDA_pooled_c <- splsda(X, y_cfactor, ncomp = 1,
                            keepX = 101)

MysPLSDA_pooled_ibd <- splsda(X, y_ibdfactor, ncomp = 1,
                              keepX = 47)

# Insepection of loading coefficients of x to see which variables were selected 
MysPLSDA_pooled_c$loadings$X
MysPLSDA_pooled_ibd$loadings$X
# The selected variables are the ones with non-zero loadings coefficients
c <- MysPLSDA_pooled_c$loadings$X[MysPLSDA_pooled_c$loadings$X !=
                                    0, ]
c
ibd <- MysPLSDA_pooled_ibd$loadings$X[MysPLSDA_pooled_ibd$loadings$X !=
                                        0, ]
ibd

# Investigate variables
cc_n <- names(c)
ibd_n <-names(ibd)

# Recode: remove repeating variables
cc_n <- gsubfn('[0-9]+', '', cc_n)
ibd_n <- gsubfn('[0-9]+', '', ibd_n)

# Remove duplicates
cc_n <- cc_n[!duplicated(cc_n)]
ibd_n <- ibd_n[!duplicated(ibd_n)]

# Intersect of variables selected by Colon cancer and IBD
int <-intersect(cc_n,ibd_n)
int

# Venndiagram
colon_lengthpls <- length(cc_n)
ibd_lengthpls <- length(ibd_n)
joint_lenthpls <-length(int)
myCol <- brewer.pal(3, "Pastel2")
myCol<-myCol[1:2]

# Plot venndiagram
pdf(file = paste0(analysis,"Vennspls.pdf"), width = 9,height = 6)
venn.plot2 <- draw.pairwise.venn(colon_lengthpls, ibd_lengthpls, joint_lenthpls, c("Colon cancer", "IBD"), 
                                 fill = myCol,
                                 lty = 'blank', 
                                 cat.pos = c(0,10), height=600, width=600, cat.cex=1, cex=1);
grid.draw(venn.plot2);
grid.newpage();
dev.off();

toc()

#----------------- 5f Random Forest -----------------#
# Original code is ukb_ML_randomforest.R by Andrea 
# This step works and takes between 3-4 hrs to run
print("Starting 5f")
tic("Step 5f")
set.seed(235)
#setwd("/rds/general/project/hda_students_data/live/Group4/General/full_scripts/data")
df <-readRDS("imputed_everything.rds")

#------------------------ Get x and y dataframes

x <- df %>% dplyr::select(-c(eid,small_intestine,rectosigmoid,rectal,anal,liver,PRS_ibd,PRS_uc,PRS_cd,cd,uc,ibd,colon))
yibd <- df %>% dplyr::select(c(ibd))
ycolon <- df %>% dplyr::select(c(colon))

#--------------------- First for IBD 100 trees

tic("IBD Random Forest Training 100 trees")
print("IBD Random Forest 100 trees")
# This takes 100-140 minutes!
rf.ibd <- randomForest(x=x, y=yibd$ibd, ntree=100, importance=TRUE)
toc()

importance<-rf.ibd$importance
#print("IBD Variable Importance as Measured by Mean Decrease Accuracy:")
#importance[order(importance[,3], decreasing =TRUE),] [1:30,]
#saveRDS(importance, "ukb_ML_rf_IBD_importance.rds")
ibd <- importance

#png("ukb_ML_rf_varImpPlotIBD.png")
#varImpPlot(rf.ibd) 
#dev.off()

#--------------------- Then for colon 100 trees

tic("Colon Random Forest Training 100 trees")
print("Colon Random Forest 100 trees")
# This takes 100-140 minutes!
rf.colon <- randomForest(x=x, y=ycolon$colon, ntree=100, importance=TRUE)
toc()

importance=rf.colon$importance
#print("Colon Variable Importance as Measured by Mean Decrease Accuracy:")
#importance[order(importance[,3], decreasing =TRUE),] [1:30,]
#saveRDS(importance, "ukb_ML_rf_Colon_importance.rds")
colon <- importance

#png("ukb_ML_rf_varImpPlotColon.png")
#varImpPlot(rf.colon) 
#dev.off()


#---------------------- Finally, the inner join

#colon <- readRDS("ukb_ML_rf_Colon_importance.rds")
#ibd <- readRDS("ukb_ML_rf_IBD_importance.rds")

# We need to pick a threshold: either the top X predictors, or a minimum MeanDecreaseAccuracy, or Gini
#colon[order(colon[,3], decreasing =TRUE),] [1:10,]
#ibd[order(ibd[,3], decreasing =TRUE),] [1:10,]

# Plots to figure out where to threshold
#library(ggplot2)
#ibd_data <- cbind(as.data.frame(ibd),rownames(ibd)) %>% transmute(ibdpredictors = rownames(ibd), MeanDecreaseAccuracy) %>% arrange(desc(MeanDecreaseAccuracy))
#ggplot(ibd_data, aes(x=reorder(ibdpredictors, MeanDecreaseAccuracy),y=MeanDecreaseAccuracy)) + geom_point()

#colon_data <- cbind(as.data.frame(colon),rownames(colon)) %>% transmute(colonpredictors = rownames(colon), MeanDecreaseAccuracy) %>% arrange(desc(MeanDecreaseAccuracy))
#ggplot(colon_data, aes(x=reorder(colonpredictors, MeanDecreaseAccuracy),y=MeanDecreaseAccuracy)) + geom_point()

# Based on the plots, there seems to be a clear drop-off, so I will set the threshold there!
threshold = 0.00025

colon_predictors <- colon %>% as.data.frame() %>% filter(MeanDecreaseAccuracy>threshold) %>% arrange(desc(MeanDecreaseAccuracy))
#print(colon_predictors)
#dim(colon_predictors)
colon_predictors <- rownames(colon_predictors)
colon_length <- length(colon_predictors)
#print("Colon Permutation Predictors")
#print(colon_predictors)
#saveRDS(colon_predictors,"ukb_ML_rf_colon_predictors.rds")

ibd_predictors <- ibd  %>% as.data.frame() %>% filter(MeanDecreaseAccuracy>threshold)%>% arrange(desc(MeanDecreaseAccuracy))
#print(ibd_predictors)
#dim(ibd_predictors)
ibd_predictors <- rownames(ibd_predictors)
ibd_length <- length(ibd_predictors)
#print("IBD Permutation Predictors")
#print(ibd_predictors)
#saveRDS(ibd_predictors,"ukb_ML_rf_IBD_predictors.rds")

#ibd_predictors_l <- c(ibd_predictors, c(0,0))
#predictors <- cbind(ibd_predictors_l,colon_predictors)
#predictors

# Then we can do an inner join on the variable names, and make a Venn diagram

print("Now we are doing the inner join")
joint_predictors <- inner_join(as.data.frame(ibd_predictors), as.data.frame(colon_predictors), by=c("ibd_predictors"="colon_predictors")) 
#print("joint predictors are: (and saving as ukb_ML_randomforest_jointpredictors_perm.rds)")
#print(joint_predictors)
saveRDS(joint_predictors,"ukb_ML_randomforest_jointpredictors_perm.rds")
#joint_length <- dim(joint_predictors)[1]
## 11 Joint predictors

toc()

#----------------- 5g LASSO -----------------#
# Original code is lasso_oem.R 
#----------------- 5h Attenuation Analysis -----------------#
# Original code is logreg_colon.R 

### These two steps are in a separate script due to very high computational demands of step 5g

print("Step 5 (part 1) complete! > just going onto Step5+ now")
toc()