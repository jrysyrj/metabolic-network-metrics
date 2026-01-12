#!/bin/bash

#SBATCH -p mit_normal
#SBATCH -o GrassmannianAngles.sh.log-%j-%a
#SBATCH -a 1-20
#SBATCH -n 5
#SBATCH --time=3:00:00

module load julia/1.9.1

MYDIR=/home/jrys/orcd/scratch/metabolic-network-metrics/organs/results
MYDIR2=/home/jrys/orcd/pool/metabolic-network-metrics/organs/results

[ ! -d "$MYDIR/organs_LNSAngles/" ] && mkdir $MYDIR/organs_LNSAngles/
[ ! -d "$MYDIR/organs_int_LNSAngles/" ] && mkdir $MYDIR/organs_int_LNSAngles/
[ ! -d "$MYDIR/organs_RNSAngles/" ] && mkdir $MYDIR/organs_RNSAngles/

#julia xx2_GrassmannianAngles.jl $SLURM_ARRAY_TASK_ID 55 $MYDIR2/organs_LNS/ $MYDIR/organs_LNSAngles/
#julia xx2_GrassmannianAngles.jl $SLURM_ARRAY_TASK_ID 55 $MYDIR2/organs_int_LNS/ $MYDIR/organs_int_LNSAngles/
julia xx2_GrassmannianAngles.jl $SLURM_ARRAY_TASK_ID 20 $MYDIR2/organs_RNS/ $MYDIR/organs_RNSAngles/
