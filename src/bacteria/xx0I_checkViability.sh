#!/bin/bash

#SBATCH -p mit_normal
#SBATCH -o checkViability.sh.log-%j
#SBATCH -n 4

module load matlab
module load gurobi

matlab -r "xx0I_checkViability_ecoli; xx0I_checkViability_ecoli_linearMOMA; xx0I_checkViability_subAGORA2; exit;"