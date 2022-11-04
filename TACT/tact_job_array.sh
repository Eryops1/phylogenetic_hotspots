#!/bin/bash
# submit_array.sh

#SBATCH --account PDiv
#SBATCH --job-name=TACT_3rd
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=melanie.tietje@bio.au.dk
#SBATCH --partition normal
#SBATCH --mem-per-cpu=150gb
#SBATCH --cpus-per-task 1
#SBATCH --time 120:00:00
#SBATCH --array=301-351


#singularity exec tact.sif tact_build_taxonomic_tree goodsp_wcp_2022_forTACT.csv --output goodsp_wcp_2022_forTACT.taxonomy.tre &&

source ~/miniconda3/bin/activate tact

echo -e "\nrunning TACT\n"


singularity exec tact.sif  tact_add_taxa --backbone gbmb_matched_monophyletic_orders.tre --taxonomy goodsp_wcp_2022_forTACT.taxonomy.tre --output output/gbmb_matched_monophyletic_orders_$SLURM_JOB_ID-$SLURM_ARRAY_TASK_ID.tacted --verbose

echo -e "\nfinished TACT\n"