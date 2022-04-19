#!/bin/bash
# submit_array.sh

#SBATCH --account PDiv
#SBATCH --job-name=nullshuffle_boss
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=melanie.tietje@bio.au.dk
#SBATCH --partition normal
#SBATCH --mem-per-cpu=1gb
#SBATCH --cpus-per-task 1
#SBATCH --time 00:05:00

rep=(1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67)
model='colwise'

export model
# write bash scripts
for ((i=0; i<=66; i++)) do
  this_rep=${rep[${i}]}
  export this_rep
  sbatch null_shuffle.sh
done
