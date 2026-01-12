#!/bin/bash

#SBATCH -p mit_normal
#SBATCH -o makeConsistent.sh.log-%j
#SBATCH -n 4

module load matlab
module load gurobi

matlab -r "xx0B_makeConsistent_ecoli_bsub;xx0B_makeConsistent;exit;"
