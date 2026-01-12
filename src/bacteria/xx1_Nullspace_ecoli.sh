#!/bin/bash

#SBATCH -p mit_normal
#SBATCH -o nullspace.sh.log-%j-%a
#SBATCH -a 1-7
#SBATCH -n 4

module load julia/1.9.1 

MYDIR=/home/jrys/orcd/pool/metabolic-network-metrics/bacteria/results
MYDIR2=/home/jrys/orcd/pool/metabolic-network-metrics/bacteria/results

[ ! -d "$MYDIR2/ecoli_LNS/" ] && mkdir $MYDIR2/ecoli_LNS/

#julia xx1_Nullspace.jl $SLURM_ARRAY_TASK_ID 7 false $MYDIR/parsed_ecoli/ $MYDIR2/ecoli_LNS/ 3156 6316

[ ! -d "$MYDIR2/ecoli_int_LNS/" ] && mkdir $MYDIR2/ecoli_int_LNS/

#julia xx1_Nullspace.jl $SLURM_ARRAY_TASK_ID 7 false $MYDIR/parsed_ecoli_int/ $MYDIR2/ecoli_int_LNS/ 3156 6316

[ ! -d "$MYDIR2/ecoli_RNS/" ] && mkdir $MYDIR2/ecoli_RNS/

julia xx1_Nullspace.jl $SLURM_ARRAY_TASK_ID 7 true $MYDIR/parsed_ecoli/ $MYDIR2/ecoli_RNS/ 3156 6316

