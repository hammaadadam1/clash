#!/bin/bash
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4G
#SBATCH --time=02:30:00
#SBATCH --array=1-1000%100
#SBATCH --output=./output/norm_%a.out

Rscript --vanilla sim_gaussian.R $SLURM_ARRAY_TASK_ID ./results
