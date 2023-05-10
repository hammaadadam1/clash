#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=01:00:00
#SBATCH --array=1-1000%200

Rscript --vanilla sim_tte.R $SLURM_ARRAY_TASK_ID ./results