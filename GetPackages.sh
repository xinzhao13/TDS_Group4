#PBS -l walltime=1:00:00
#PBS -l select=1:ncpus=1:mem=1gb

path=/rds/general/project/hda_students_data/live/Group4/General/full_scripts
cd $path

module load anaconda3/personal
module load fix_unwritable_tmp
date

# STEP 5 : 

Rscript GetPackages.R


date