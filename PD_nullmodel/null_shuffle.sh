#!/bin/bash

#SBATCH --account PDiv
#SBATCH --job-name=nullshuffle
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=melanie.tietje@bio.au.dk
#SBATCH --partition normal
#SBATCH --mem-per-cpu=8gb
#SBATCH --cpus-per-task 1
#SBATCH --time 00:60:00
#SBATCH --output=/dev/null

source ~/miniconda3/bin/activate R-env-4
Rscript null_shuffle.R $this_rep $model >logfile.txt
