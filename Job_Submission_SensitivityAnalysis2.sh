#!/bin/bash
#PBS -l walltime=30:00:00
#PBS -l select=1:ncpus=16:mem=10gb
#PBS -N gwas_imp_yayy
#PBS -J 1-22

# Please set the directory of the R Scripts here, which will be carried throughout the entire code:
path=/rds/general/project/hda_students_data/live/Group4/General/full_scripts
cd $path

module load plink

geno_path=/rds/general/project/uk-biobank-2018/live/reference/sdata_12032018
fam_path=$path/data
data_path=$path/data
results_path=$path/data/v4
snps_path=$path/data

plink --bfile $geno_path/ukb_imp_chr$PBS_ARRAY_INDEX \
--fam $fam_path/0-2-ukb_imp_yayy.fam \
--extract $snps_path/0-2-ibd_snps_list.txt \
--covar $data_path/0-1-confounders.txt keep-pheno-on-missing-cov --no-const-covar \
--maf 0.001 \
--ci 0.95 \
--covar-number 1-12 \
--hide-covar \
--logistic beta \
--out $results_path/logistic_imp_v4_chr$PBS_ARRAY_INDEX