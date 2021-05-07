#----------------- 6.2: Sensitivity Analysis & Assumptions Check -----------------#

# This script performs the actions of Sensitivity Analysis & Assumptions Check
# It is designed to run on the HPC server

print("initiating Step 6.2 now (Sensitivity Analysis)")

# It requires
# > all previous steps
# > UKB genotype data
# > Several packages
library(tidyverse)
library(tictoc)
library(data.table)
library(MendelianRandomization)
library(patchwork)
# Documentation: https://cran.r-project.org/web/packages/MendelianRandomization/MendelianRandomization.pdf
# Github repo: https://github.com/cran/MendelianRandomization

tic("Step 6.2")

args=commandArgs(trailingOnly=TRUE)
nchunks=as.numeric(args[1])
path=as.character(args[2])
print(path)
datapath<-paste0(path,"/data/")
analysis<-paste0(path,"/analysis/")
setwd(datapath)
print("We've set the working directory")

# It outputs analysis and visualisation in the folder "analysis"

#----------------- 6d Run mini GWAS -----------------#
# Original code is 0-3-logistic_geno.sh
##### This has just been completed
print("Step 6d complete")

#----------------- 6e Consolidate mini GWAS output from 22 chromosomes -----------------#
# Original code is 0-4-consolidate_gwas.R
print("Starting 6e")
tic("Step 6e")

# Read in our gwas regressions
files<-list.files(path="v4", pattern="*.assoc.logistic", full.names=TRUE, recursive=FALSE)
cc_gwas_summary<-do.call("rbind",lapply(files, fread))

write.csv(cc_gwas_summary,"0-4-cc_gwas_v4.csv",row.names=FALSE)

toc()

#----------------- 6f Produce all methods sensitivity analysis -----------------#
# Original code is 1-1-mr_sensitivity.R
print("Starting 6f")
tic("Step 6f")

# Define paths and set working directory
y_path<-datapath
X_path<-datapath
Z_path<-datapath

# Set working directory
setwd(y_path)

### Data prep -----------------------------------------------------

# Read in data
X<-readRDS(paste0(X_path,"corrected_betas_ibd.rds"))  # beta=logOR and corrected, no se
y<-read.csv(paste0(y_path,"0-4-cc_gwas_v4.csv")) # beta=logOR
Z<-read.csv(paste0(Z_path,"iv_maf_ibd.txt")) # beta=logOR, has se

# Read in se from Z for X that's not in corrected_betas_ibd.rds, but required by mr_obj() function
X<-left_join(X,Z[,c("rsid","se")],by=c("id"="rsid"))

X<-X[order(X$id),]
y<-y[order(y$SNP),]

# Review top rows
#head(X)
#head(y)

# Extract list of SNPs for later
X_snps<-X$id

# Create MR object
mr_obj<-mr_input(bx=X$newbeta,
                 bxse=X$se,
                 by=y$BETA,
                 byse=y$SE)

# Sensitivity analysis --------------------------------------------------

# All methods => stats
mr_all<-mr_allmethods(mr_obj,method="all")

# Write to file
sink(file="1-1-mr_all_stats.txt",type="output",append=TRUE)
print(mr_all)
sink() # stop writing to file



## mendelian randomisation mr_plot() function has issue today 11/04/2021
## Below are my attempt to reconstruct the plot -----------------------------

df<-mr_all@Values
n<-nrow(df)

df_new<-data.frame(matrix(NA, nrow=n, ncol=3))

df_new<-df[c(2,4,8), 1:2]
df_new$Intercept<-c(rep(0,2), df[9, 2])
df_new$Method<-c("WM","IVW","MR-Egger")


betaX=mr_all@Data@betaX
betaY=mr_all@Data@betaY
exposure="IBD"
outcome="colon cancer"

# Plot the SNPs
plt_scatter<-
  ggplot(data = NULL, aes(x=betaX, y=betaY)) +
  geom_point(color="grey34",size=1) +
  geom_hline(yintercept = 0, color = "grey23", alpha = 0.2) +
  geom_vline(xintercept = 0, color = "grey23", alpha = 0.2) +
  
  # Plot lines of standard MR + sensitivity tests
  geom_abline(data = df_new, aes(intercept = Intercept, slope = Estimate, 
                                 color = Method, linetype = Method),
              show.legend = TRUE, size = 1) +
  
  scale_colour_manual(name="Method", labels = df_new$Method,
                      breaks = c("WM",  "IVW", "MR-Egger"),
                      values = c("#F8766D", "#69b3a2", "#404080")) +
  scale_linetype_manual(name="Method", labels = df_new$Method,
                        breaks = c("WM",  "IVW", "MR-Egger"),
                        values = c("solid", "solid", "solid")) +
  xlab(paste("Genetic association with", exposure)) +
  ylab(paste("Genetic association with", outcome)) + theme_classic()


# Forest plot -----------------------------------------------------------


# Extract from mr_obj object
bx<-mr_obj@betaX
by<-mr_obj@betaY
bxse<-mr_obj@betaXse
byse<-mr_obj@betaYse
snps<-mr_obj@snps

alpha<-0.05

estimates<-by/bx
ci_range<-qnorm(1-alpha/2)

ci_lower<-estimates-(ci_range*byse)/abs(bx)
ci_upper<-estimates+(ci_range*byse)/abs(bx)

#Create the dataframe with these values
df<-data.frame(snps, estimates, ci_lower, ci_upper, color=rep(0,78), names=rep(0,78))

egger_output<-mr_egger(mr_obj, alpha=alpha) # <- actually we are plotting MR-Egger estimates now
egger_estimate<-egger_output$Estimate
egger_ci_lower<-egger_output$CILower.Est
egger_ci_upper<-egger_output$CIUpper.Est

ivw_output<-mr_ivw(mr_obj, alpha=alpha) # <- actually we are plotting MR-ivw estimates now
ivw_estimate<-ivw_output$Estimate
ivw_ci_lower<-ivw_output$CILower
ivw_ci_upper<-ivw_output$CIUpper

wm_output<-mr_median(mr_obj, alpha=alpha) # <- actually we are plotting MR-wm estimates now
wm_estimate<-wm_output$Estimate
wm_ci_lower<-wm_output$CILower
wm_ci_upper<-wm_output$CIUpper

# Derive egger and ivw outputs - default used in mr_forest() from MendelianRandomisation package
egger_row<-data.frame("MR-Egger estimate", egger_estimate, egger_ci_lower, egger_ci_upper,color="tomato", `names`="MR-Egger")
names(egger_row)<-names(df)

ivw_row<-data.frame("IVW estimate", ivw_estimate, ivw_ci_lower, ivw_ci_upper,color="tomato", `names`="IVW")
names(ivw_row)<-names(df)

wm_row<-data.frame("wm estimate", wm_estimate, wm_ci_lower, wm_ci_upper,color="tomato", `names`="WM")
names(wm_row)<-names(df)


df$color<-ifelse(df$ci_lower*df$ci_upper<0,"grey39","blue1")
df$names<-as.character(X$id)

#The idea is that we first order the estimates, then reverse the order
#Then add the calculated ones to the beginning which is actually the end
factor_order<-rev(df$snps[order(df$estimates)])
df$snps<-factor(df$snps, levels = factor_order)

factor_order_names<-rev(df$names[order(df$estimates)])
df$names<-factor(df$names, levels = factor_order_names)

# x<-c(ivw_row$snps, factor_order)
df<-rbind(df, wm_row, ivw_row, egger_row)


# Original function from the package has issue here as input was already factored
df$snps<-factor(df$snps, levels= c(as.character(egger_row$snps), 
                                   as.character(ivw_row$snps), as.character(wm_row$snps),as.character(factor_order)))

x_label <- "IBD-associated genetic variants"
interval_type <- paste(100*(1-alpha), "% CI)", sep = "")
y_label <-  paste("Causal estimate on colon cancer (", interval_type, sep = "")

# Forest: appendix --------------------------------------------------------

# Export to png 
# png("1-1-mr_plot_forest.jpg", width = 11.69, height = 8.27, units = 'in', res=600)
plt_forest_appendix<-ggplot(data = df, aes(y=snps, x=estimates, xmin=ci_lower, xmax = ci_upper)) +
  geom_point(shape = 15) +
  geom_point(data = ivw_row, aes(y=snps, x = estimates), shape = 18, size = 3, color="tomato") +
  geom_linerange(color=df$color) +
  geom_vline(xintercept = 0, lty= 2) +
  ylab(x_label) + xlab(y_label) +
  theme_update(axis.text.y = element_text(color = "grey20",
                                          size = 5, angle = 0, hjust = 1, vjust = 0))+
  scale_y_discrete(labels=c("MR-Egger Estimate", "IVW Estimate", "WM Estimate",levels(df$names)))+
  theme_classic()

# Forest: results --------------------------------------------------------

# Export to png 
# png("1-1-mr_plot_forest.jpg", width = 11.69, height = 8.27, units = 'in', res=600)
plt_forest_results<-ggplot(data = df, aes(y=snps, x=estimates, xmin=ci_lower, xmax = ci_upper)) +
  geom_point(shape = 15) +
  geom_point(data = ivw_row, aes(y=snps, x = estimates), shape = 18, size = 3, color="tomato") +
  geom_linerange(color=df$color) +
  geom_vline(xintercept = 0, lty= 2) +
  ylab(x_label) + xlab(y_label) +
  theme_update(axis.text.y = element_text(color = "grey20",
                                          size = 5, angle = 0, hjust = 1, vjust = 0))+
  scale_y_discrete(labels=c(rep("",81)))+
  theme_classic() +
  theme(axis.text.y=element_text(size=5))

### Figure exports -----------------------------------------------------------

patchwork<-plt_forest_results / plt_scatter
patch_results<-patchwork + plot_annotation(tag_levels = list(c('a)','b)')))


# Export to png
png(paste0(analysis,"results; forest+scatter.jpg"), width = 8.69, height = 8.27, units = 'in', res=600)
patch_results
dev.off()

# Export to png
png(paste0(analysis,"appendix; forest.jpg"), width = 6, height = 8.27, units = 'in', res=600)
plt_forest_appendix
dev.off()

# Export to png
png(paste0(analysis,"appendix; scatter.jpg"), width = 11, height = 8.27, units = 'in', res=600)
plt_scatter
dev.off()


toc()

#----------------- 6g Generate PhenoScanner output -----------------#
# Original code is 1-2-phenoscanner.R
print("Starting 6g")
tic("Step 6g")

# Define paths and set working directory
y_path<-datapath
X_path<-datapath

### Data prep -----------------------------------------------------

# Read in data
X<-readRDS("corrected_betas_ibd.rds")  # beta=logOR and corrected, no se

# Extract list of SNPs for later
X_snps<-X$id
# Pheno scanner -----------------------------------------------------

# Export phenoscanner output to file
for (i in (1:length(X_snps))){
  snp <- X_snps[i]
  snp_data<-phenoscanner(snp)
  write.table(snp_data,paste0(analysis,"1-2-phenoscanner_output.csv"),
              append=TRUE,sep="/")
  
  print(paste0("Progress: ",i,"/",length(X_snps)))  
}

toc()


print("Step 6 fully complete!")
toc()