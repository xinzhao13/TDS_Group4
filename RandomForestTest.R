#----------------- 5: Analysis & Visualisation -----------------#

# This script is complete and takes about 28+RF mins to run

# This script performs the actions of Analysis & Visualisation, defined by this flowchart: https://whimsical.com/tds-r-scripts-and-data-flow-VmAm6BzY1jUML2a32569t2
# It is designed to run on the HPC server

# It requires
# > ukb_hes_everything.rds
# > Several packages
library(tidyverse)#installed
library(tictoc)#installed
library(gtsummary)#installed
library(gt)#installed
library(data.table)#installed
library(patchwork)#installed
library(ggplot2)#installed
library(plotROC)#installed
library(rcompanion) #for nagelkerke #installed
library(imputeMissings)#installed
library(sgPLS)#installed
library(VennDiagram)#installed
library(pheatmap)#installed
library(utils)#installed
library(ROSE) #installed
library(gsubfn) #installed
library(randomForest)#installed
library(gridExtra) #installed
library(Hmisc)#installed
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

#----------------- 5f Random Forest -----------------#
# Original code is ukb_ML_randomforest.R by Andrea 
# This step works and takes between XX & XX mins to run
print("Starting 5f")
tic("Step 5f")

#setwd("/rds/general/project/hda_students_data/live/Group4/General/full_scripts/data")
df <-readRDS("imputed_everything.rds")

#------------------------ Get x and y dataframes

x <- df %>% dplyr::select(-c(eid,small_intestine,rectosigmoid,rectal,anal,liver,PRS_ibd,PRS_uc,PRS_cd,cd,uc,ibd,colon))
yibd <- df %>% dplyr::select(c(ibd))
ycolon <- df %>% dplyr::select(c(colon))

#--------------------- First for IBD 100 trees

tic("IBD Random Forest Training 100 trees")
print("IBD Random Forest 100 trees")
# 1 tree takes 83 seconds, so 100 trees will take 2hrs 20 mins
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
# 1 tree takes 83 seconds, so 100 trees will take 2hrs 20 mins
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

# I picked an arbitrary value of 230 as the Gini threshold, as this should work fine no matter the balance between case and control:

# Plots to figure out where to threshold
#library(ggplot2)
#ibd_data <- cbind(as.data.frame(ibd),rownames(ibd)) %>% transmute(ibdpredictors = rownames(ibd), MeanDecreaseAccuracy) %>% arrange(desc(MeanDecreaseAccuracy))
#ggplot(ibd_data, aes(x=reorder(ibdpredictors, MeanDecreaseAccuracy),y=MeanDecreaseAccuracy)) + geom_point()

#colon_data <- cbind(as.data.frame(colon),rownames(colon)) %>% transmute(colonpredictors = rownames(colon), MeanDecreaseAccuracy) %>% arrange(desc(MeanDecreaseAccuracy))
#ggplot(colon_data, aes(x=reorder(colonpredictors, MeanDecreaseAccuracy),y=MeanDecreaseAccuracy)) + geom_point()

# Based on the plots, there seems to be a clear drop-off, so I will set the threshold there!
threshold = 0.0003

colon_gini_predictors <- colon %>% as.data.frame() %>% filter(MeanDecreaseAccuracy>threshold) %>% arrange(desc(MeanDecreaseAccuracy))
#print(colon_gini_predictors)
#dim(colon_gini_predictors)
colon_gini_predictors <- rownames(colon_gini_predictors)
colon_length <- length(colon_gini_predictors)

ibd_gini_predictors <- ibd  %>% as.data.frame() %>% filter(MeanDecreaseAccuracy>threshold)%>% arrange(desc(MeanDecreaseAccuracy))
#print(ibd_gini_predictors)
#dim(ibd_gini_predictors)
ibd_gini_predictors <- rownames(ibd_gini_predictors)
ibd_length <- length(ibd_gini_predictors)

#ibd_gini_predictors_l <- c(ibd_gini_predictors, c(0,0))
#predictors <- cbind(ibd_gini_predictors_l,colon_gini_predictors)
#predictors

# Then we can do an inner join on the variable names, and make a Venn diagram

print("Now we are doing the inner join")
joint_predictors <- inner_join(as.data.frame(ibd_gini_predictors), as.data.frame(colon_gini_predictors), by=c("ibd_gini_predictors"="colon_gini_predictors")) 
#print("joint predictors are: (and saving as ukb_ML_randomforest_jointpredictors_perm.rds)")
#print(joint_predictors)
saveRDS(joint_predictors,"ukb_ML_randomforest_jointpredictors_perm.rds")
joint_length <- dim(joint_predictors)[1]
## 11 Joint predictors

toc()

#----------------- 5g LASSO -----------------#
# Original code is lasso_oem.R by Marie 
#----------------- 5h Attenuation Analysis -----------------#
# Original code is logreg_colon.R by Marie 

### These two steps are in a separate script due to very high computational demands of step 5g

print("Step 5 (part 1) complete! > just going onto Step5+ now")
toc()