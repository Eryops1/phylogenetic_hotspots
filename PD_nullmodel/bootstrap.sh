#!/bin/bash

#SBATCH --account PDiv
#SBATCH --job-name=subsample
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=melanie.tietje@bio.au.dk
#SBATCH --partition normal
<<<<<<< HEAD
#SBATCH --mem-per-cpu=64gb
=======
#SBATCH --mem-per-cpu=32gb
>>>>>>> 0571433986b7bd0045f2427c8893bd3144bda85a
#SBATCH --cpus-per-task 1
#SBATCH --time 00:20:00
##SBATCH --output=/dev/null

source ~/miniconda3/bin/activate R-env-4
Rscript manual_bootstrap.R $this_rep > bootlog.txt
