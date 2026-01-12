#!/bin/bash

#SBATCH -p mit_normal
#SBATCH -o GrassmannianAngles.sh.log-%j-%a
#SBATCH -a 1-15
#SBATCH -n 2

module load julia/1.9.1

MYDIR=/home/jrys/orcd/scratch/metabolic-network-metrics/bacteria/results
MYDIR2=/home/jrys/orcd/pool/metabolic-network-metrics/bacteria/results

[ ! -d "$MYDIR/ecoli_LNSAngles/" ] && mkdir $MYDIR/ecoli_LNSAngles/
[ ! -d "$MYDIR/ecoli_int_LNSAngles/" ] && mkdir $MYDIR/ecoli_int_LNSAngles/
[ ! -d "$MYDIR/ecoli_RNSAngles/" ] && mkdir $MYDIR/ecoli_RNSAngles/

#julia xx2_GrassmannianAngles.jl $SLURM_ARRAY_TASK_ID 15 $MYDIR2/ecoli_LNS/ $MYDIR/ecoli_LNSAngles/

#julia xx2_GrassmannianAngles.jl $SLURM_ARRAY_TASK_ID 15 $MYDIR2/ecoli_int_LNS/ $MYDIR/ecoli_int_LNSAngles/

julia xx2_GrassmannianAngles.jl $SLURM_ARRAY_TASK_ID 15 $MYDIR2/ecoli_RNS/ $MYDIR/ecoli_RNSAngles/
