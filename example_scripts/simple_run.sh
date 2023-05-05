#!/bin/bash
#SBATCH -t 3-00:00
#SBATCH --job-name= ER
#SBATCH --mail-user=
#SBATCH --mail-type=FAIL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=150g
#SBATCH --cpus-per-task=16

module load gcc/10.2.0
module load r/4.2.0

Rscript simple_run.R --yaml_path 'path_to_yaml' --coarse_grid T
