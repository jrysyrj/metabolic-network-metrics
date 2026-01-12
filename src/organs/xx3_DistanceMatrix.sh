#!/bin/bash

#SBATCH -p mit_normal
#SBATCH -o DistanceMatrix.sh.log-%j-%a
#SBATCH -a 1-20
#SBATCH -n 5
#SBATCH --time=5:00:00

module load julia/1.9.1

MYDIR=/home/jrys/orcd/pool/metabolic-network-metrics/organs/results
MYDIR2=/home/jrys/orcd/scratch/metabolic-network-metrics/organs/results

[ ! -d "$MYDIR/distances/" ] && mkdir $MYDIR/distances/

#julia xx3_DistanceMatrix.jl $SLURM_ARRAY_TASK_ID 55 $MYDIR2/organs_LNSAngles/ $MYDIR/organs_LNS/ $MYDIR/distances/organs_LNS
#julia xx3_DistanceMatrix.jl $SLURM_ARRAY_TASK_ID 55 $MYDIR2/organs_int_LNSAngles/ $MYDIR/organs_int_LNS/ $MYDIR/distances/organs_int_LNS
julia xx3_DistanceMatrix.jl $SLURM_ARRAY_TASK_ID 20 $MYDIR2/organs_RNSAngles/ $MYDIR/organs_RNS/ $MYDIR/distances/organs_RNS

