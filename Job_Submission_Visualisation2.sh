#PBS -l walltime=8:00:00
#PBS -l select=1:ncpus=10:mem=120gb

# Please set the directory of the R Scripts here, which will be carried throughout the entire code:
path=/rds/general/project/hda_students_data/live/Group4/General/full_scripts
cd $path

nchunks=10

module load anaconda3/personal
module load fix_unwritable_tmp

date

# STEP 5plus :

Rscript 5plus_Analysis_and_Visualisation.R $nchunks $path

date