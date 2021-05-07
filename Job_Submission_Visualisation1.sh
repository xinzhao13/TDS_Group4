#PBS -l walltime=6:00:00
#PBS -l select=1:ncpus=1:mem=20gb

# This script takes approx 4 hours to run from start to finish

# Please set the directory of the R Scripts here, which will be carried throughout each R script:
path=/rds/general/project/hda_students_data/live/Group4/General/full_scripts
cd $path

nchunks=1

module load anaconda3/personal
module load fix_unwritable_tmp
date

# STEP 5 :

Rscript 5_Analysis_and_Visualisation.R $nchunks $path

date