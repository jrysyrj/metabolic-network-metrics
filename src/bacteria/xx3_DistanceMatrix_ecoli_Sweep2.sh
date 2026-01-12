#!/bin/bash

#SBATCH -p mit_normal
#SBATCH -o DistanceMatrix.sh.log-%j-%a
#SBATCH -a 1-55
#SBATCH -n 1
#SBATCH --time=6:00:00

module load julia/1.9.1

MYDIR=/home/jrys/orcd/pool/metabolic-network-metrics/bacteria/results
MYDIR2=/home/jrys/orcd/scratch/metabolic-network-metrics/bacteria/results

[ ! -d "$MYDIR/distances/" ] && mkdir $MYDIR/distances/
[ ! -d "$MYDIR/distances/ecoli_nullspace_sweep2/" ] && mkdir $MYDIR/distances/ecoli_nullspace_sweep2/


julia xx3_DistanceMatrix.jl $SLURM_ARRAY_TASK_ID 55 $MYDIR2/ecoli_int_LNS2Angles_1/ $MYDIR2/ecoli_int_LNS2/1/ $MYDIR/distances/ecoli_nullspace_sweep2/ecoli_int_LNS2_1
julia xx3_DistanceMatrix.jl $SLURM_ARRAY_TASK_ID 55 $MYDIR2/ecoli_int_LNS2Angles_2/ $MYDIR2/ecoli_int_LNS2/2/ $MYDIR/distances/ecoli_nullspace_sweep2/ecoli_int_LNS2_2
julia xx3_DistanceMatrix.jl $SLURM_ARRAY_TASK_ID 55 $MYDIR2/ecoli_int_LNS2Angles_3/ $MYDIR2/ecoli_int_LNS2/3/ $MYDIR/distances/ecoli_nullspace_sweep2/ecoli_int_LNS2_3
julia xx3_DistanceMatrix.jl $SLURM_ARRAY_TASK_ID 55 $MYDIR2/ecoli_int_LNS2Angles_4/ $MYDIR2/ecoli_int_LNS2/4/ $MYDIR/distances/ecoli_nullspace_sweep2/ecoli_int_LNS2_4

#julia xx3_DistanceMatrix.jl $SLURM_ARRAY_TASK_ID 55 $MYDIR2/ecoli_LNS2Angles_1/ $MYDIR2/ecoli_LNS2/1/ $MYDIR/distances/ecoli_nullspace_sweep2/ecoli_LNS2_1
#julia xx3_DistanceMatrix.jl $SLURM_ARRAY_TASK_ID 55 $MYDIR2/ecoli_LNS2Angles_2/ $MYDIR2/ecoli_LNS2/2/ $MYDIR/distances/ecoli_nullspace_sweep2/ecoli_LNS2_2
#julia xx3_DistanceMatrix.jl $SLURM_ARRAY_TASK_ID 55 $MYDIR2/ecoli_LNS2Angles_3/ $MYDIR2/ecoli_LNS2/3/ $MYDIR/distances/ecoli_nullspace_sweep2/ecoli_LNS2_3
#julia xx3_DistanceMatrix.jl $SLURM_ARRAY_TASK_ID 55 $MYDIR2/ecoli_LNS2Angles_4/ $MYDIR2/ecoli_LNS2/4/ $MYDIR/distances/ecoli_nullspace_sweep2/ecoli_LNS2_4

#julia xx3_DistanceMatrix.jl $SLURM_ARRAY_TASK_ID 55 $MYDIR2/ecoli_RNS2Angles_1/ $MYDIR2/ecoli_RNS2/1/ $MYDIR/distances/ecoli_nullspace_sweep2/ecoli_RNS2_1
#julia xx3_DistanceMatrix.jl $SLURM_ARRAY_TASK_ID 55 $MYDIR2/ecoli_RNS2Angles_2/ $MYDIR2/ecoli_RNS2/2/ $MYDIR/distances/ecoli_nullspace_sweep2/ecoli_RNS2_2
#julia xx3_DistanceMatrix.jl $SLURM_ARRAY_TASK_ID 55 $MYDIR2/ecoli_RNS2Angles_3/ $MYDIR2/ecoli_RNS2/3/ $MYDIR/distances/ecoli_nullspace_sweep2/ecoli_RNS2_3
#julia xx3_DistanceMatrix.jl $SLURM_ARRAY_TASK_ID 55 $MYDIR2/ecoli_RNS2Angles_4/ $MYDIR2/ecoli_RNS2/4/ $MYDIR/distances/ecoli_nullspace_sweep2/ecoli_RNS2_4
