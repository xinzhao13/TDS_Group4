#PBS -l walltime=6:00:00
#PBS -l select=1:ncpus=24:mem=30gb
#PBS -N 1nodesevcores

cd /rds/general/project/hda_students_data/live/Group4/General/full_scripts/
module load anaconda3/personal
module load fix_unwritable_tmp

nchunks=24
Rscript 1_SNP_List_Creation.R $nchunks