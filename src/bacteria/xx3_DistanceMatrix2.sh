#!/bin/bash

#SBATCH -p mit_normal
#SBATCH -o DistanceMatrix2.sh.log-%j-%a
#SBATCH -a 1-45
#SBATCH -n 2

module load julia/1.9.1

MYDIR=/home/jrys/orcd/pool/metabolic-network-metrics/bacteria/results
MYDIR2=/home/jrys/orcd/scratch/metabolic-network-metrics/bacteria/results

[ ! -d "$MYDIR/distances/" ] && mkdir $MYDIR/distances/

#julia xx3_DistanceMatrix2.jl $SLURM_ARRAY_TASK_ID 45 $MYDIR2/ecoli_subAGORA2_LNSAngles/ $MYDIR/ecoli_LNS/ $MYDIR/subAGORA2_LNS/ $MYDIR/distances/ecoli_subAGORA2_LNS

julia xx3_DistanceMatrix2.jl $SLURM_ARRAY_TASK_ID 45 $MYDIR2/ecoli_subAGORA2_int_LNSAngles/ $MYDIR/ecoli_int_LNS/ $MYDIR/subAGORA2_int_LNS/ $MYDIR/distances/ecoli_subAGORA2_int_LNS

#julia xx3_DistanceMatrix2.jl $SLURM_ARRAY_TASK_ID 45 $MYDIR2/ecoli_subAGORA2_RNSAngles/ $MYDIR/ecoli_RNS/ $MYDIR/subAGORA2_RNS/ $MYDIR/distances/ecoli_subAGORA2_RNS

