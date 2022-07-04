#!/bin/bash

#SBATCH --account PDiv
#SBATCH --job-name=ses_mpd
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=melanie.tietje@bio.au.dk
#SBATCH --partition normal
#SBATCH --mem-per-cpu=16gb
#SBATCH --cpus-per-task 1
#SBATCH --time 2:00:00
#SBATCH --output=ses_mpd.out

source ~/miniconda3/bin/activate R-env-4
Rscript ses_mpd.R >ses_mpd.txt
