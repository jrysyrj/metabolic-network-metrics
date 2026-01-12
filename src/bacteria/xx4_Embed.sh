#!/bin/bash

#SBATCH -o embed.sh.log-%j-%a
#SBATCH -a 1-20
#SBATCH -n 5

module load julia/1.9.1

MYDIR=/home/jrys/orcd/pool/metabolic-network-metrics/bacteria/results/

[ ! -d "$MYDIR/sgdMMDS/" ] && mkdir $MYDIR/sgdMMDS/
[ ! -d "$MYDIR/sgdMMDS/LGRASSMANN_subAGORA2" ] && mkdir $MYDIR/sgdMMDS/LGRASSMANN_subAGORA2
[ ! -d "$MYDIR/sgdMMDS/LGRASSMANN_subAGORA2_int" ] && mkdir $MYDIR/sgdMMDS/LGRASSMANN_subAGORA2_int
[ ! -d "$MYDIR/sgdMMDS/RGRASSMANN_subAGORA2" ] && mkdir $MYDIR/sgdMMDS/RGRASSMANN_subAGORA2

[ ! -d "$MYDIR/sgdMMDS/LGRASSMANN_subAGORA2/stress" ] && mkdir $MYDIR/sgdMMDS/LGRASSMANN_subAGORA2/stress
[ ! -d "$MYDIR/sgdMMDS/LGRASSMANN_subAGORA2_int/stress" ] && mkdir $MYDIR/sgdMMDS/LGRASSMANN_subAGORA2_int/stress
[ ! -d "$MYDIR/sgdMMDS/RGRASSMANN_subAGORA2/stress" ] && mkdir $MYDIR/sgdMMDS/RGRASSMANN_subAGORA2/stress

[ ! -d "$MYDIR/sgdMMDS/LGRASSMANN_subAGORA2/distortion" ] && mkdir $MYDIR/sgdMMDS/LGRASSMANN_subAGORA2/distortion
[ ! -d "$MYDIR/sgdMMDS/LGRASSMANN_subAGORA2_int/distortion" ] && mkdir $MYDIR/sgdMMDS/LGRASSMANN_subAGORA2_int/distortion
[ ! -d "$MYDIR/sgdMMDS/RGRASSMANN_subAGORA2/distortion" ] && mkdir $MYDIR/sgdMMDS/RGRASSMANN_subAGORA2/distortion

julia xx4_Embed.jl $SLURM_ARRAY_TASK_ID 20
