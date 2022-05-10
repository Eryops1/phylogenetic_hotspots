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

vars=('hfp' 'deforest' 'bio1' 'bio5' 'bio6' 'bio7' 'bio12' 'bio15' 'PC_primf' 'PC_primn' 'PC_secdf' 'PC_secdn' 'PC_urban' 'PC_c3ann' 'PC_c4ann' 'PC_c3per' 'PC_c4per' 'PC_c3nfx' 'PC_pastr' 'PC_range' 'PC_secmb' 'PC_secma' 'FC_primf' 'FC_primn' 'FC_secdf' 'FC_secdn' 'FC_urban' 'FC_c3ann' 'FC_c4ann' 'FC_c3per' 'FC_c4per' 'FC_c3nfx' 'FC_pastr' 'FC_range' 'FC_secmb' 'FC_secma')

# pass on variables to child scripts
for ((i=0; i<=35; i++)) do
  this_var=${vars[${i}]}
  export this_var
  sbatch env.sh
done
