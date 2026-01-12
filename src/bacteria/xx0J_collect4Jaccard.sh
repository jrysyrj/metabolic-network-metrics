#!/bin/bash

#SBATCH -p mit_normal
#SBATCH -o collect4Jaccard.sh.log-%j
#SBATCH -n 4

module load matlab
module load gurobi

matlab -r "xx0J_collect4Jaccard; exit;"