#!/bin/bash

#SBATCH -p mit_normal
#SBATCH -o nullspace.sh.log-%j-%a
#SBATCH -a 1-20
#SBATCH -n 2

module load julia/1.9.1 

MYDIR=/home/jrys/orcd/pool/metabolic-network-metrics/bacteria/results
MYDIR2=/home/jrys/orcd/scratch/metabolic-network-metrics/bacteria/results


[ ! -d "$MYDIR2/ecoli_LNS/" ] && mkdir $MYDIR2/ecoli_LNS/ && mkdir $MYDIR2/ecoli_LNS/1/ && mkdir $MYDIR2/ecoli_LNS/2/ && mkdir $MYDIR2/ecoli_LNS/3/ && mkdir $MYDIR2/ecoli_LNS/4/
#[ ! -d "$MYDIR2/ecoli_RNS/" ] && mkdir $MYDIR2/ecoli_RNS/ && mkdir $MYDIR2/ecoli_RNS/1/ && mkdir $MYDIR2/ecoli_RNS/2/ && mkdir $MYDIR2/ecoli_RNS/3/ && mkdir $MYDIR2/ecoli_RNS/4/

julia xx1_NullspaceSweep.jl $SLURM_ARRAY_TASK_ID 20 false $MYDIR/parsed_ecoli/ $MYDIR2/ecoli_LNS/ 3156 6316
#julia xx1_NullspaceSweep.jl $SLURM_ARRAY_TASK_ID 20 true $MYDIR/parsed_ecoli/ $MYDIR2/ecoli_RNS/ 3156 6316

