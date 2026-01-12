#!/bin/bash

#SBATCH -p mit_normal
#SBATCH -o GrassmannianAngles.sh.log-%j-%a
#SBATCH -a 1-55
#SBATCH -n 1
#SBATCH --time=6:00:00

module load julia/1.9.1

MYDIR=/home/jrys/orcd/scratch/metabolic-network-metrics/bacteria/results

[ ! -d "$MYDIR/ecoli_int_LNS2Angles_1/" ] && mkdir $MYDIR/ecoli_int_LNS2Angles_1/
#[ ! -d "$MYDIR/ecoli_LNS2Angles_1/" ] && mkdir $MYDIR/ecoli_LNS2Angles_1/
#[ ! -d "$MYDIR/ecoli_RNS2Angles_1/" ] && mkdir $MYDIR/ecoli_RNS2Angles_1/

[ ! -d "$MYDIR/ecoli_int_LNS2Angles_2/" ] && mkdir $MYDIR/ecoli_int_LNS2Angles_2/
#[ ! -d "$MYDIR/ecoli_LNS2Angles_2/" ] && mkdir $MYDIR/ecoli_LNS2Angles_2/
#[ ! -d "$MYDIR/ecoli_RNS2Angles_2/" ] && mkdir $MYDIR/ecoli_RNS2Angles_2/

[ ! -d "$MYDIR/ecoli_int_LNS2Angles_3/" ] && mkdir $MYDIR/ecoli_int_LNS2Angles_3/
#[ ! -d "$MYDIR/ecoli_LNS2Angles_3/" ] && mkdir $MYDIR/ecoli_LNS2Angles_3/
#[ ! -d "$MYDIR/ecoli_RNS2Angles_3/" ] && mkdir $MYDIR/ecoli_RNS2Angles_3/

[ ! -d "$MYDIR/ecoli_int_LNS2Angles_4/" ] && mkdir $MYDIR/ecoli_int_LNS2Angles_4/
#[ ! -d "$MYDIR/ecoli_LNS2Angles_4/" ] && mkdir $MYDIR/ecoli_LNS2Angles_4/
#[ ! -d "$MYDIR/ecoli_RNS2Angles_4/" ] && mkdir $MYDIR/ecoli_RNS2Angles_4/

#julia xx2_GrassmannianAngles.jl $SLURM_ARRAY_TASK_ID 55 $MYDIR/ecoli_int_LNS2/1/ $MYDIR/ecoli_int_LNS2Angles_1/
#julia xx2_GrassmannianAngles.jl $SLURM_ARRAY_TASK_ID 55 $MYDIR/ecoli_LNS2/1/ $MYDIR/ecoli_LNS2Angles_1/
#julia xx2_GrassmannianAngles.jl $SLURM_ARRAY_TASK_ID 55 $MYDIR/ecoli_RNS2/1/ $MYDIR/ecoli_RNS2Angles_1/

#julia xx2_GrassmannianAngles.jl $SLURM_ARRAY_TASK_ID 55 $MYDIR/ecoli_int_LNS2/2/ $MYDIR/ecoli_int_LNS2Angles_2/
#julia xx2_GrassmannianAngles.jl $SLURM_ARRAY_TASK_ID 55 $MYDIR/ecoli_LNS2/2/ $MYDIR/ecoli_LNS2Angles_2/
#julia xx2_GrassmannianAngles.jl $SLURM_ARRAY_TASK_ID 55 $MYDIR/ecoli_RNS2/2/ $MYDIR/ecoli_RNS2Angles_2/

julia xx2_GrassmannianAngles.jl $SLURM_ARRAY_TASK_ID 55 $MYDIR/ecoli_int_LNS2/3/ $MYDIR/ecoli_int_LNS2Angles_3/
#julia xx2_GrassmannianAngles.jl $SLURM_ARRAY_TASK_ID 55 $MYDIR/ecoli_LNS2/3/ $MYDIR/ecoli_LNS2Angles_3/
#julia xx2_GrassmannianAngles.jl $SLURM_ARRAY_TASK_ID 55 $MYDIR/ecoli_RNS2/3/ $MYDIR/ecoli_RNS2Angles_3/

julia xx2_GrassmannianAngles.jl $SLURM_ARRAY_TASK_ID 55 $MYDIR/ecoli_int_LNS2/4/ $MYDIR/ecoli_int_LNS2Angles_4/
#julia xx2_GrassmannianAngles.jl $SLURM_ARRAY_TASK_ID 55 $MYDIR/ecoli_LNS2/4/ $MYDIR/ecoli_LNS2Angles_4/
#julia xx2_GrassmannianAngles.jl $SLURM_ARRAY_TASK_ID 55 $MYDIR/ecoli_RNS2/4/ $MYDIR/ecoli_RNS2Angles_4/
