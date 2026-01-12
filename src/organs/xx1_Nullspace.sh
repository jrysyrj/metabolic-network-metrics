#!/bin/bash

#SBATCH -p mit_normal
#SBATCH -o nullspace.sh.log-%j-%a
#SBATCH -a 1-20
#SBATCH -n 5
#SBATCH --time=2:00:00

module load julia/1.9.1 

MYDIR=/home/jrys/orcd/pool/metabolic-network-metrics/organs/results
MYDIR2=/home/jrys/orcd/pool/metabolic-network-metrics/organs/results

#[ ! -d "$MYDIR2/organs_LNS/" ] && mkdir $MYDIR2/organs_LNS/
#julia xx1_Nullspace.jl $SLURM_ARRAY_TASK_ID 5 false $MYDIR/parsed_organs/ $MYDIR2/organs_LNS/ 10882 20561
#[ ! -d "$MYDIR2/organs_int_LNS/" ] && mkdir $MYDIR2/organs_int_LNS/
#julia xx1_Nullspace.jl $SLURM_ARRAY_TASK_ID 5 false $MYDIR/parsed_organs_int/ $MYDIR2/organs_int_LNS/ 10882 20561
[ ! -d "$MYDIR2/organs_RNS/" ] && mkdir $MYDIR2/organs_RNS/
julia xx1_Nullspace.jl $SLURM_ARRAY_TASK_ID 20 true $MYDIR/parsed_organs/ $MYDIR2/organs_RNS/ 10882 20561

