#!/bin/bash
# submit_array.sh

#SBATCH --account PDiv
#SBATCH --job-name=tipshuffle_master
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=melanie.tietje@bio.au.dk
#SBATCH --partition normal
#SBATCH --mem-per-cpu=1gb
#SBATCH --cpus-per-task 1
#SBATCH --time 00:05:00

rep=seq 67

export rep
# write bash scripts
for ((i=0; i<=66; i++)) do
  this_rep=${rep[${i}]}
  export this_rep
  sbatch tipshuffle.sh
done
