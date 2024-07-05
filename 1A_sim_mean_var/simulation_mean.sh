#!/bin/bash 
#SBATCH --job-name=sim_mean
#SBATCH --array 1-300
#SBATCH --output out/change_%A_%a.out
#SBATCH --error out/change_%A_%a.err
#SBATCH --time=24:00:00
#SBATCH --partition=medium
#SBATCH --ntasks=1
#SBATCH --mem=4G
#SBATCH --mail-type=FAIL

module purge
module load r/4.3.0

R CMD BATCH --vanilla "simulation_mean.R" "out/$SLURM_ARRAY_JOB_ID-$SLURM_ARRAY_TASK_ID.Routput"
