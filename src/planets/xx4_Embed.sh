#!/bin/bash

#SBATCH -o embed.sh.log-%j-%a
#SBATCH -a 1-20
#SBATCH -n 5
#SBATCH --time=5:00:00

module load julia/1.9.1

MYDIR=/home/jrys/orcd/pool/metabolic-network-metrics/planets/results/

[ ! -d "$MYDIR/sgdMMDS/" ] && mkdir $MYDIR/sgdMMDS/
[ ! -d "$MYDIR/sgdMMDS/LGRASSMANN_planets" ] && mkdir $MYDIR/sgdMMDS/LGRASSMANN_planets
[ ! -d "$MYDIR/sgdMMDS/LGRASSMANN_planets_int" ] && mkdir $MYDIR/sgdMMDS/LGRASSMANN_planets_int
[ ! -d "$MYDIR/sgdMMDS/RGRASSMANN_planets" ] && mkdir $MYDIR/sgdMMDS/RGRASSMANN_planets

[ ! -d "$MYDIR/sgdMMDS/LGRASSMANN_planets/stress" ] && mkdir $MYDIR/sgdMMDS/LGRASSMANN_planets/stress
[ ! -d "$MYDIR/sgdMMDS/LGRASSMANN_planets_int/stress" ] && mkdir $MYDIR/sgdMMDS/LGRASSMANN_planets_int/stress
[ ! -d "$MYDIR/sgdMMDS/RGRASSMANN_planets/stress" ] && mkdir $MYDIR/sgdMMDS/RGRASSMANN_planets/stress

[ ! -d "$MYDIR/sgdMMDS/LGRASSMANN_planets/distortion" ] && mkdir $MYDIR/sgdMMDS/LGRASSMANN_planets/distortion
[ ! -d "$MYDIR/sgdMMDS/LGRASSMANN_planets_int/distortion" ] && mkdir $MYDIR/sgdMMDS/LGRASSMANN_planets_int/distortion
[ ! -d "$MYDIR/sgdMMDS/RGRASSMANN_planets/distortion" ] && mkdir $MYDIR/sgdMMDS/RGRASSMANN_planets/distortion

julia xx4_Embed.jl $SLURM_ARRAY_TASK_ID 20
