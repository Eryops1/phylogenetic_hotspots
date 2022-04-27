#!/bin/bash
# submit_array.sh

#SBATCH --account PDiv
#SBATCH --job-name=env_var_boss
#SBATCH --mail-type=FAIL,END
##SBATCH --mail-user=melanie.tietje@bio.au.dk
#SBATCH --partition normal
#SBATCH --mem-per-cpu=1gb
#SBATCH --cpus-per-task 1
#SBATCH --time 00:05:00

vars=('hfp', 'deforest', 'CHELSA_bio1', 'CHELSA_bio5', 'CHELSA_bio6', 'CHELSA_bio7', 'CHELSA_bio12')

# pass on variables to child scripts
for ((i=0; i<=6; i++)) do
  this_var=${vars[${i}]}
  export this_var
  sbatch env.sh
done
