#!/bin/bash

#SBATCH -p mit_normal
#SBATCH -o GrassmannianAngles.sh.log-%j-%a
#SBATCH -a 1-55
#SBATCH -n 1

module load julia/1.9.1

MYDIR=/home/jrys/orcd/scratch/metabolic-network-metrics/bacteria/results

[ ! -d "$MYDIR/ecoli_LNSAngles_1/" ] && mkdir $MYDIR/ecoli_LNSAngles_1/
#[ ! -d "$MYDIR/ecoli_RNSAngles_1/" ] && mkdir $MYDIR/ecoli_RNSAngles_1/

[ ! -d "$MYDIR/ecoli_LNSAngles_2/" ] && mkdir $MYDIR/ecoli_LNSAngles_2/
#[ ! -d "$MYDIR/ecoli_RNSAngles_2/" ] && mkdir $MYDIR/ecoli_RNSAngles_2/

[ ! -d "$MYDIR/ecoli_LNSAngles_3/" ] && mkdir $MYDIR/ecoli_LNSAngles_3/
#[ ! -d "$MYDIR/ecoli_RNSAngles_3/" ] && mkdir $MYDIR/ecoli_RNSAngles_3/

[ ! -d "$MYDIR/ecoli_LNSAngles_4/" ] && mkdir $MYDIR/ecoli_LNSAngles_4/
#[ ! -d "$MYDIR/ecoli_RNSAngles_4/" ] && mkdir $MYDIR/ecoli_RNSAngles_4/

julia xx2_GrassmannianAngles.jl $SLURM_ARRAY_TASK_ID 55 $MYDIR/ecoli_LNS/1/ $MYDIR/ecoli_LNSAngles_1/
#julia xx2_GrassmannianAngles.jl $SLURM_ARRAY_TASK_ID 55 $MYDIR/ecoli_RNS/1/ $MYDIR/ecoli_RNSAngles_1/

julia xx2_GrassmannianAngles.jl $SLURM_ARRAY_TASK_ID 55 $MYDIR/ecoli_LNS/2/ $MYDIR/ecoli_LNSAngles_2/
#julia xx2_GrassmannianAngles.jl $SLURM_ARRAY_TASK_ID 55 $MYDIR/ecoli_RNS/2/ $MYDIR/ecoli_RNSAngles_2/

julia xx2_GrassmannianAngles.jl $SLURM_ARRAY_TASK_ID 55 $MYDIR/ecoli_LNS/3/ $MYDIR/ecoli_LNSAngles_3/
#julia xx2_GrassmannianAngles.jl $SLURM_ARRAY_TASK_ID 55 $MYDIR/ecoli_RNS/3/ $MYDIR/ecoli_RNSAngles_3/

julia xx2_GrassmannianAngles.jl $SLURM_ARRAY_TASK_ID 55 $MYDIR/ecoli_LNS/4/ $MYDIR/ecoli_LNSAngles_4/
#julia xx2_GrassmannianAngles.jl $SLURM_ARRAY_TASK_ID 55 $MYDIR/ecoli_RNS/4/ $MYDIR/ecoli_RNSAngles_4/
