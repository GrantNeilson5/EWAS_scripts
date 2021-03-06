#!/bin/sh
#PBS -V #export all enviroment variables to the batch job
#PBS -q pq #submit to mrc high memory Queue, this can be changed on ISCA wiki  
#PBS -l walltime=24:00:00 # maximum wall time for the job
#PBS -A Research_Project-MRC190311 #reserch project to submit under
#PBS -l feature=highmem
#PBS -l nodes=1:ppn=16 # Number of processors 
#PBS -m e -M g.neilson@exeter.ac.uk # email me at job completion 

module load Miniconda3

source activate QC_packages  # activating environemt 

cd /gpfs/mrc0/projects/Research_Project-MRC193462/EWAS/EWAS_scripts/ # Setting Working directory to where script is located 

Rscript EWAS_analysis_in_parralel.R ## execute R script
