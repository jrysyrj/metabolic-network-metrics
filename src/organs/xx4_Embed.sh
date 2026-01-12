#!/bin/bash

#SBATCH -o embed.sh.log-%j-%a
#SBATCH -a 1-20
#SBATCH -n 5
#SBATCH --time=5:00:00

module load julia/1.9.1

MYDIR=/home/jrys/orcd/pool/metabolic-network-metrics/organs/results/

[ ! -d "$MYDIR/sgdMMDS/" ] && mkdir $MYDIR/sgdMMDS/
[ ! -d "$MYDIR/sgdMMDS/LGRASSMANN_organs" ] && mkdir $MYDIR/sgdMMDS/LGRASSMANN_organs
[ ! -d "$MYDIR/sgdMMDS/LGRASSMANN_organs_int" ] && mkdir $MYDIR/sgdMMDS/LGRASSMANN_organs_int
[ ! -d "$MYDIR/sgdMMDS/RGRASSMANN_organs" ] && mkdir $MYDIR/sgdMMDS/RGRASSMANN_organs

[ ! -d "$MYDIR/sgdMMDS/LGRASSMANN_organs/stress" ] && mkdir $MYDIR/sgdMMDS/LGRASSMANN_organs/stress
[ ! -d "$MYDIR/sgdMMDS/LGRASSMANN_organs_int/stress" ] && mkdir $MYDIR/sgdMMDS/LGRASSMANN_organs_int/stress
[ ! -d "$MYDIR/sgdMMDS/RGRASSMANN_organs/stress" ] && mkdir $MYDIR/sgdMMDS/RGRASSMANN_organs/stress

[ ! -d "$MYDIR/sgdMMDS/LGRASSMANN_organs/distortion" ] && mkdir $MYDIR/sgdMMDS/LGRASSMANN_organs/distortion
[ ! -d "$MYDIR/sgdMMDS/LGRASSMANN_organs_int/distortion" ] && mkdir $MYDIR/sgdMMDS/LGRASSMANN_organs_int/distortion
[ ! -d "$MYDIR/sgdMMDS/RGRASSMANN_organs/distortion" ] && mkdir $MYDIR/sgdMMDS/RGRASSMANN_organs/distortion

julia xx4_Embed.jl $SLURM_ARRAY_TASK_ID 20
