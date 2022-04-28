#!/bin/bash

#SBATCH --account PDiv
#SBATCH --job-name=env_var_child
#SBATCH --mail-type=FAIL,END
##SBATCH --mail-user=melanie.tietje@bio.au.dk
#SBATCH --partition normal
#SBATCH --mem-per-cpu=32gb
#SBATCH --cpus-per-task 1
#SBATCH --time 02:00:00
##SBATCH --output=NULL

source ~/miniconda3/bin/activate pdiv
Rscript env_var.R $this_var > logfile_"$thisvar".txt
