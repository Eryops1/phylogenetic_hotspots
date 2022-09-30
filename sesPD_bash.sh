#!/bin/bash

#SBATCH --account PDiv
#SBATCH --job-name=picante
#SBATCH --job-name=picante_SES_PD
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=melanie.tietje@bio.au.dk
#SBATCH --partition normal
#SBATCH --mem-per-cpu=50gb
#SBATCH --cpus-per-task 1
#SBATCH --time 168:00:00
#SBATCH --output=picante.out

source ~/miniconda3/bin/activate R-env-4

Rscript sesPD_picante.R

