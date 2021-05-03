#PBS -l walltime=1:00:00
#PBS -l select=1:ncpus=1:mem=10gb

# Please set the directory of the R Scripts here, which will be carried throughout the entire code:
path=/rds/general/project/hda_students_data/live/Group4/General/full_scripts
cd $path

module load anaconda3/personal
module load fix_unwritable_tmp
date

# STEP 4 : 

Rscript 4_Final_Data_Creation.R $path

date