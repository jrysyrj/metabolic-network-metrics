#!/bin/bash

#SBATCH -p mit_normal
#SBATCH -o nullspace.sh.log-%j-%a
#SBATCH -a 1-5
#SBATCH -n 1
#SBATCH -N 1

module load julia/1.9.1 

MYDIR=/home/jrys/orcd/pool/metabolic-network-metrics/bacteria/results
MYDIR2=/home/jrys/orcd/pool/metabolic-network-metrics/bacteria/results

[ ! -d "$MYDIR2/subAGORA2_LNS/" ] && mkdir $MYDIR2/subAGORA2_LNS/

julia xx1_Nullspace.jl $SLURM_ARRAY_TASK_ID 5 false $MYDIR/parsed_subAGORA2/ $MYDIR2/subAGORA2_LNS/ 3156 6316

[ ! -d "$MYDIR2/subAGORA2_int_LNS/" ] && mkdir $MYDIR2/subAGORA2_int_LNS/

julia xx1_Nullspace.jl $SLURM_ARRAY_TASK_ID 5 false $MYDIR/parsed_subAGORA2_int/ $MYDIR2/subAGORA2_int_LNS/ 3156 6316

[ ! -d "$MYDIR2/subAGORA2_RNS/" ] && mkdir $MYDIR2/subAGORA2_RNS/

julia xx1_Nullspace.jl $SLURM_ARRAY_TASK_ID 5 true $MYDIR/parsed_subAGORA2/ $MYDIR2/subAGORA2_RNS/ 3156 6316

