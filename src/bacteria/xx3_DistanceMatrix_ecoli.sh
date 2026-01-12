#!/bin/bash

#SBATCH -p mit_normal
#SBATCH -o DistanceMatrix.sh.log-%j-%a
#SBATCH -a 1-15
#SBATCH -n 1

module load julia/1.9.1

MYDIR=/home/jrys/orcd/pool/metabolic-network-metrics/bacteria/results
MYDIR2=/home/jrys/orcd/scratch/metabolic-network-metrics/bacteria/results

[ ! -d "$MYDIR/distances/" ] && mkdir $MYDIR/distances/

#julia xx3_DistanceMatrix.jl $SLURM_ARRAY_TASK_ID 15 $MYDIR2/ecoli_LNSAngles/ $MYDIR/ecoli_LNS/ $MYDIR/distances/ecoli_LNS
#julia xx3_DistanceMatrix.jl $SLURM_ARRAY_TASK_ID 15 $MYDIR2/ecoli_int_LNSAngles/ $MYDIR/ecoli_int_LNS/ $MYDIR/distances/ecoli_int_LNS
julia xx3_DistanceMatrix.jl $SLURM_ARRAY_TASK_ID 15 $MYDIR2/ecoli_RNSAngles/ $MYDIR/ecoli_RNS/ $MYDIR/distances/ecoli_RNS

