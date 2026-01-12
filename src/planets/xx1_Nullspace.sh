#!/bin/bash

#SBATCH -p mit_normal
#SBATCH -o nullspace.sh.log-%j-%a
#SBATCH -a 1-5
#SBATCH -n 5
#SBATCH --time=2:00:00

module load julia/1.9.1 

MYDIR=/home/jrys/orcd/pool/metabolic-network-metrics/planets/results
MYDIR2=/home/jrys/orcd/pool/metabolic-network-metrics/planets/results

[ ! -d "$MYDIR2/planets_LNS/" ] && mkdir $MYDIR2/planets_LNS/
julia xx1_Nullspace.jl $SLURM_ARRAY_TASK_ID 5 false $MYDIR/parsed_planets/ $MYDIR2/planets_LNS/ 216 722
[ ! -d "$MYDIR2/planets_int_LNS/" ] && mkdir $MYDIR2/planets_int_LNS/
julia xx1_Nullspace.jl $SLURM_ARRAY_TASK_ID 5 false $MYDIR/parsed_planets_int/ $MYDIR2/planets_int_LNS/ 216 722
[ ! -d "$MYDIR2/planets_RNS/" ] && mkdir $MYDIR2/planets_RNS/
julia xx1_Nullspace.jl $SLURM_ARRAY_TASK_ID 5 true $MYDIR/parsed_planets/ $MYDIR2/planets_RNS/ 216 722

