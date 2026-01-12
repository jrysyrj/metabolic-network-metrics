#!/bin/bash

#SBATCH -p mit_normal
#SBATCH -o GrassmannianAngles2.sh.log-%j-%a
#SBATCH -a 1-40
#SBATCH -n 1

module load julia/1.9.1

MYDIR=/home/jrys/orcd/scratch/metabolic-network-metrics/bacteria/results
MYDIR2=/home/jrys/orcd/pool/metabolic-network-metrics/bacteria/results

[ ! -d "$MYDIR/ecoli_subAGORA2_LNSAngles/" ] && mkdir $MYDIR/ecoli_subAGORA2_LNSAngles/
[ ! -d "$MYDIR/ecoli_subAGORA2_int_LNSAngles/" ] && mkdir $MYDIR/ecoli_subAGORA2_int_LNSAngles/
[ ! -d "$MYDIR/ecoli_subAGORA2_RNSAngles/" ] && mkdir $MYDIR/ecoli_subAGORA2_RNSAngles/

julia xx2_GrassmannianAngles2.jl $SLURM_ARRAY_TASK_ID 40 $MYDIR2/ecoli_LNS/ $MYDIR2/subAGORA2_LNS/ $MYDIR/ecoli_subAGORA2_LNSAngles/ 

julia xx2_GrassmannianAngles2.jl $SLURM_ARRAY_TASK_ID 40 $MYDIR2/ecoli_int_LNS/ $MYDIR2/subAGORA2_int_LNS/ $MYDIR/ecoli_subAGORA2_int_LNSAngles/

julia xx2_GrassmannianAngles2.jl $SLURM_ARRAY_TASK_ID 40 $MYDIR2/ecoli_RNS/ $MYDIR2/subAGORA2_RNS/ $MYDIR/ecoli_subAGORA2_RNSAngles/