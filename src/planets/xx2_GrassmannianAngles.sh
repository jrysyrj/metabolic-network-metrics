#!/bin/bash

#SBATCH -p mit_normal
#SBATCH -o GrassmannianAngles.sh.log-%j-%a
#SBATCH -a 1-10
#SBATCH -n 5
#SBATCH --time=3:00:00

module load julia/1.9.1

MYDIR=/home/jrys/orcd/scratch/metabolic-network-metrics/planets/results
MYDIR2=/home/jrys/orcd/pool/metabolic-network-metrics/planets/results

[ ! -d "$MYDIR/planets_LNSAngles/" ] && mkdir $MYDIR/planets_LNSAngles/
[ ! -d "$MYDIR/planets_int_LNSAngles/" ] && mkdir $MYDIR/planets_int_LNSAngles/
[ ! -d "$MYDIR/planets_RNSAngles/" ] && mkdir $MYDIR/planets_RNSAngles/

julia xx2_GrassmannianAngles.jl $SLURM_ARRAY_TASK_ID 10 $MYDIR2/planets_LNS/ $MYDIR/planets_LNSAngles/
julia xx2_GrassmannianAngles.jl $SLURM_ARRAY_TASK_ID 10 $MYDIR2/planets_int_LNS/ $MYDIR/planets_int_LNSAngles/
julia xx2_GrassmannianAngles.jl $SLURM_ARRAY_TASK_ID 10 $MYDIR2/planets_RNS/ $MYDIR/planets_RNSAngles/
