#PBS -l walltime=5:00:00
#PBS -l select=1:ncpus=24:mem=80gb
#PBS -N 1nodesevcores

cd /rds/general/project/hda_students_data/live/Group4/General/full_scripts/
module load anaconda3/personal
module load fix_unwritable_tmp

# STEP 2: Genetic Data Creation - approx 7 mins runtime

date

# import_snps.sh
module load plink
cd /rds/general/project/hda_students_data/live/Group4/General/full_scripts/data/import_snps_results_all/

UKB=/rds/general/project/uk-biobank-2018/live/reference/sdata_12032018
FAM=/rds/general/project/hda_students_data/live/Group4/General/full_scripts/data
LIST=/rds/general/project/hda_students_data/live/Group4/General/full_scripts/data
for chr in $(seq 1 22)
do
echo $chr

plink --bfile $UKB/ukb_imp_chr$chr --fam $FAM/ukb_imp.fam --extract $LIST/snps_imp_list_all.txt --make-bed --out ukb_imp_selected_chr$chr

done

# merge_snps.sh
module load plink
cd /rds/general/project/hda_students_data/live/Group4/General/full_scripts/data/

CHRL=/rds/general/project/hda_students_data/live/Group4/General/full_scripts/data
RES=/rds/general/project/hda_students_data/live/Group4/General/full_scripts/data/import_snps_results_all

for chr in $(cat $CHRL/chr_list.txt)
do cat $RES/ukb_imp_selected_chr${chr}.bim; done > $RES/ukb_imp_merged.bim

(echo -en "\x6C\x1B\x01"; for chr in $(cat $CHRL/chr_list.txt)
do tail -c +4 $RES/ukb_imp_selected_chr${chr}.bed; done) > $RES/ukb_imp_merged.bed

cd /rds/general/project/hda_students_data/live/Group4/General/full_scripts/

date

nchunks=24
Rscript 2_Genetic_Data_Creation.R $nchunks

date