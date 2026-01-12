#!/bin/bash

#SBATCH -p mit_normal
#SBATCH -o DistanceMatrix.sh.log-%j-%a
#SBATCH -a 1-55
#SBATCH -n 1

module load julia/1.9.1

MYDIR=/home/jrys/orcd/pool/metabolic-network-metrics/bacteria/results
MYDIR2=/home/jrys/orcd/scratch/metabolic-network-metrics/bacteria/results

[ ! -d "$MYDIR/distances/" ] && mkdir $MYDIR/distances/
[ ! -d "$MYDIR/distances/ecoli_nullspace_sweep/" ] && mkdir $MYDIR/distances/ecoli_nullspace_sweep/

julia xx3_DistanceMatrix.jl $SLURM_ARRAY_TASK_ID 55 $MYDIR2/ecoli_LNSAngles_1/ $MYDIR2/ecoli_LNS/1/ $MYDIR/distances/ecoli_nullspace_sweep/ecoli_LNS_1
julia xx3_DistanceMatrix.jl $SLURM_ARRAY_TASK_ID 55 $MYDIR2/ecoli_LNSAngles_2/ $MYDIR2/ecoli_LNS/2/ $MYDIR/distances/ecoli_nullspace_sweep/ecoli_LNS_2
julia xx3_DistanceMatrix.jl $SLURM_ARRAY_TASK_ID 55 $MYDIR2/ecoli_LNSAngles_3/ $MYDIR2/ecoli_LNS/3/ $MYDIR/distances/ecoli_nullspace_sweep/ecoli_LNS_3
julia xx3_DistanceMatrix.jl $SLURM_ARRAY_TASK_ID 55 $MYDIR2/ecoli_LNSAngles_4/ $MYDIR2/ecoli_LNS/4/ $MYDIR/distances/ecoli_nullspace_sweep/ecoli_LNS_4

#julia xx3_DistanceMatrix.jl $SLURM_ARRAY_TASK_ID 55 $MYDIR2/ecoli_RNSAngles_1/ $MYDIR2/ecoli_RNS/1/ $MYDIR/distances/ecoli_nullspace_sweep/ecoli_RNS_1
#julia xx3_DistanceMatrix.jl $SLURM_ARRAY_TASK_ID 55 $MYDIR2/ecoli_RNSAngles_2/ $MYDIR2/ecoli_RNS/2/ $MYDIR/distances/ecoli_nullspace_sweep/ecoli_RNS_2
#julia xx3_DistanceMatrix.jl $SLURM_ARRAY_TASK_ID 55 $MYDIR2/ecoli_RNSAngles_3/ $MYDIR2/ecoli_RNS/3/ $MYDIR/distances/ecoli_nullspace_sweep/ecoli_RNS_3
#julia xx3_DistanceMatrix.jl $SLURM_ARRAY_TASK_ID 55 $MYDIR2/ecoli_RNSAngles_4/ $MYDIR2/ecoli_RNS/4/ $MYDIR/distances/ecoli_nullspace_sweep/ecoli_RNS_4
