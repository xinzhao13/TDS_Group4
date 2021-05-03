#PBS -l walltime=1:00:00
#PBS -l select=1:ncpus=1:mem=30gb
#PBS -N 1nodesevcores

cd /rds/general/project/hda_students_data/live/Group4/General/full_scripts/
module load anaconda3/personal
module load fix_unwritable_tmp
date

# STEP 3

Rscript 3_Outcome_Covariate_Data_Creation.R

date