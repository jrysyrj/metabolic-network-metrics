#!/bin/bash

#SBATCH -p mit_normal
#SBATCH -o ecoliKOs.sh.log-%j
#SBATCH -n 4

module load matlab
module load gurobi

matlab -r "xx0G_ecoliKOs;exit;"