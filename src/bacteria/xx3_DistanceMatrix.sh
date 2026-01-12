#!/bin/bash

#SBATCH -p mit_normal
#SBATCH -o DistanceMatrix.sh.log-%j-%a
#SBATCH -a 1-15
#SBATCH -n 1

module load julia/1.9.1

MYDIR=/home/jrys/orcd/pool/metabolic-network-metrics/bacteria/results
MYDIR2=/home/jrys/orcd/scratch/metabolic-network-metrics/bacteria/results

[ ! -d "$MYDIR/distances/" ] && mkdir $MYDIR/distances/

#julia xx3_DistanceMatrix.jl $SLURM_ARRAY_TASK_ID 15 $MYDIR2/subAGORA2_LNSAngles/ $MYDIR/subAGORA2_LNS/ $MYDIR/distances/subAGORA2_LNS
#julia xx3_DistanceMatrix.jl $SLURM_ARRAY_TASK_ID 15 $MYDIR2/subAGORA2_int_LNSAngles/ $MYDIR/subAGORA2_int_LNS/ $MYDIR/distances/subAGORA2_int_LNS
julia xx3_DistanceMatrix.jl $SLURM_ARRAY_TASK_ID 15 $MYDIR2/subAGORA2_RNSAngles/ $MYDIR/subAGORA2_RNS/ $MYDIR/distances/subAGORA2_RNS

