#!/bin/bash

#SBATCH -p mit_normal
#SBATCH -o DistanceMatrix.sh.log-%j-%a
#SBATCH -a 1-10
#SBATCH -n 5
#SBATCH --time=5:00:00

module load julia/1.9.1

MYDIR=/home/jrys/orcd/pool/metabolic-network-metrics/planets/results
MYDIR2=/home/jrys/orcd/scratch/metabolic-network-metrics/planets/results

[ ! -d "$MYDIR/distances/" ] && mkdir $MYDIR/distances/

julia xx3_DistanceMatrix.jl $SLURM_ARRAY_TASK_ID 10 $MYDIR2/planets_LNSAngles/ $MYDIR/planets_LNS/ $MYDIR/distances/planets_LNS
julia xx3_DistanceMatrix.jl $SLURM_ARRAY_TASK_ID 10 $MYDIR2/planets_int_LNSAngles/ $MYDIR/planets_int_LNS/ $MYDIR/distances/planets_int_LNS
julia xx3_DistanceMatrix.jl $SLURM_ARRAY_TASK_ID 10 $MYDIR2/planets_RNSAngles/ $MYDIR/planets_RNS/ $MYDIR/distances/planets_RNS

