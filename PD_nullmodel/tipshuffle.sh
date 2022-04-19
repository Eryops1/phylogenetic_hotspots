#!/bin/bash

#SBATCH --account PDiv
#SBATCH --job-name=tipshuffle
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=melanie.tietje@bio.au.dk
#SBATCH --partition normal
#SBATCH --mem-per-cpu=8gb
#SBATCH --cpus-per-task 1
#SBATCH --time 00:30:00
#SBATCH --output=/dev/null

source ~/miniconda3/bin/activate R-env-4
Rscript tipshuffle.R $this_rep > log_"$SLURM_JOB_ID"_"$this_rep".txt
