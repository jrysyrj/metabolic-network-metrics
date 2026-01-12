#!/bin/bash

#SBATCH -p mit_normal
#SBATCH -o GrassmannianAngles.sh.log-%j-%a
#SBATCH -a 1-20
#SBATCH -n 1

module load julia/1.9.1

MYDIR=/home/jrys/orcd/scratch/metabolic-network-metrics/bacteria/results
MYDIR2=/home/jrys/orcd/pool/metabolic-network-metrics/bacteria/results

#[ ! -d "$MYDIR/subAGORA2_LNSAngles/" ] && mkdir $MYDIR/subAGORA2_LNSAngles/
#[ ! -d "$MYDIR/subAGORA2_int_LNSAngles/" ] && mkdir $MYDIR/subAGORA2_int_LNSAngles/
[ ! -d "$MYDIR/subAGORA2_RNSAngles/" ] && mkdir $MYDIR/subAGORA2_RNSAngles/

#julia xx2_GrassmannianAngles.jl $SLURM_ARRAY_TASK_ID 20 $MYDIR2/subAGORA2_LNS/ $MYDIR/subAGORA2_LNSAngles/
#julia xx2_GrassmannianAngles.jl $SLURM_ARRAY_TASK_ID 20 $MYDIR2/subAGORA2_int_LNS/ $MYDIR/subAGORA2_int_LNSAngles/

julia xx2_GrassmannianAngles.jl $SLURM_ARRAY_TASK_ID 20 $MYDIR2/subAGORA2_RNS/ $MYDIR/subAGORA2_RNSAngles/
