#PBS -l walltime=4:00:00
#PBS -l select=1:ncpus=12:mem=120gb

# Please set the directory of the R Scripts here, which will be carried throughout the entire code:
path=/rds/general/project/hda_students_data/live/Group4/General/full_scripts
cd $path

nchunks=12

module load anaconda3/personal
module load fix_unwritable_tmp
date

# STEP : 

Rscript 5plus_Analysis_and_Visualisation.R $nchunks $path

date