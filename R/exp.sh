#!/bin/bash

#SBATCH --job-name=array
#SBATCH --output=EfficientCCA/logs/array_missing_%A_%a.out
#SBATCH --error=EfficientCCA/logs/array_missing_%A_%a.err
#SBATCH --array=1-30
#SBATCH --time=35:00:00
#SBATCH --partition=caslake
#SBATCH --cpus-per-task=8
#SBATCH --ntasks=1
#SBATCH --mem=5G
#SBATCH --account=pi-cdonnat


echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "My SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID
# Add lines here to run your computations
job_id=$SLURM_ARRAY_JOB_ID
#module load libgmp
module load R/4.2.0

result_file="new_exp_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
echo "result file is ${result_file}"

Rscript EfficientCCA/R/run_experiments.R $SLURM_ARRAY_TASK_ID $result_file $1 $2 $3
