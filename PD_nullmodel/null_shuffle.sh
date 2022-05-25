#!/bin/bash

#SBATCH --account PDiv
#SBATCH --job-name=all_indices
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=melanie.tietje@bio.au.dk
#SBATCH --partition normal
#SBATCH --mem-per-cpu=64gb
#SBATCH --cpus-per-task 1
#SBATCH --time 7:00:00
#SBATCH --output=nullshuffle_indices.out

source ~/miniconda3/bin/activate R-env-4
Rscript null_shuffle.R $this_rep $this_model >logfile_14h.txt
