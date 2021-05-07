#PBS -l walltime=6:00:00
#PBS -l select=1:ncpus=1:mem=20gb

# Please set the directory of the R Scripts here, which will be carried throughout each R script:
path=/rds/general/project/hda_students_data/live/Group4/General/full_scripts
cd $path

nchunks=1

# It takes less than one hour to run.

module load anaconda3/personal
module load fix_unwritable_tmp
date

# STEP 6 : Sensitivity Analysis

Rscript 6.1_Sensitivity_Analysis.R $nchunks $path

date