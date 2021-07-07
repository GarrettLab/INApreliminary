#!/bin/bash
#SBATCH --account=epi
#SBATCH --qos=epi
#SBATCH --job-name=INA  #Job name
#SBATCH --mail-type=ALL   # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=<username>@ufl.edu  # Where to send mail
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32gb   # Per processor memory
#SBATCH --time=01:00:00   # Walltime
#SBATCH --output=ina_%j.out   # Name output file
# Displya date host and directory
date; hostname; pwd

# loading modules
module load ufrc
module load R

#load files
Rscript --vanilla INA_exe.R
