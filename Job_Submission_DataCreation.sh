#PBS -l walltime=3:00:00
#PBS -l select=1:ncpus=24:mem=80gb

# Please set the directory of the R Scripts here, which will be carried throughout each R script:
path=/rds/general/project/hda_students_data/live/Group4/General/full_scripts
cd $path

# It takes 80-120 mins to run from start to finish.

module load anaconda3/personal
module load fix_unwritable_tmp
date

# STEP 1 : SNP List Creation - approx 60-90 mins runtime

nchunks=24
Rscript 1_SNP_List_Creation.R $nchunks $path

date

# STEP 2: Genetic Data Creation - approx 7 mins runtime

# import_snps.sh
module load plink
cd $path/data/import_snps_results_all/

UKB=/rds/general/project/uk-biobank-2018/live/reference/sdata_12032018
FAM=$path/data
LIST=$path/data
for chr in $(seq 1 22)
do
echo $chr

plink --bfile $UKB/ukb_imp_chr$chr --fam $FAM/ukb_imp.fam --extract $LIST/snps_imp_list_all.txt --make-bed --out ukb_imp_selected_chr$chr

done

# merge_snps.sh
module load plink
cd $path/data/

CHRL=$path/data
RES=$path/data/import_snps_results_all

for chr in $(cat $CHRL/chr_list.txt)
do cat $RES/ukb_imp_selected_chr${chr}.bim; done > $RES/ukb_imp_merged.bim

(echo -en "\x6C\x1B\x01"; for chr in $(cat $CHRL/chr_list.txt)
do tail -c +4 $RES/ukb_imp_selected_chr${chr}.bed; done) > $RES/ukb_imp_merged.bed

cd $path

date

nchunks=24
Rscript 2_Genetic_Data_Creation.R $nchunks $path

date

# STEP 3 : Outcome & Covariate Data Creation - approx 25 mins runtime

Rscript 3_Outcome_Covariate_Data_Creation.R $path

date

# STEP 4 : Final Data Creation - approx 10 mins runtime

Rscript 4_Final_Data_Creation.R $path

date