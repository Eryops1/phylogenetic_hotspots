#!/bin/bash

#SBATCH --account PDiv
<<<<<<< HEAD
#SBATCH --job-name=picante
=======
#SBATCH --job-name=picante SES_PD
>>>>>>> eb3e0b88f1189ea00a5abcc0d9462d2527741ad8
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=melanie.tietje@bio.au.dk
#SBATCH --partition normal
#SBATCH --mem-per-cpu=50gb
#SBATCH --cpus-per-task 1
#SBATCH --time 168:00:00
#SBATCH --output=picante.out

source ~/miniconda3/bin/activate R-env-4

Rscript sesPD_picante.R

