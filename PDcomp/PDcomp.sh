#!/bin/bash

#SBATCH --account PDiv
#SBATCH --job-name=PDcomp
##SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=melanie.tietje@bio.au.dk
#SBATCH --partition normal
#SBATCH --mem-per-cpu=8gb
#SBATCH --cpus-per-task 1
#SBATCH --time 00:15:00
#SBATCH --output=NULL

source ~/miniconda3/bin/activate R-env-4
Rscript 01c_get_PDcomplementarity_cluster.R $this_rep
