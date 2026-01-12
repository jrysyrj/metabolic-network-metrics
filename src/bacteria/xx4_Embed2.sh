#!/bin/bash

#SBATCH -o embed2.sh.log-%j-%a
#SBATCH -a 1-20
#SBATCH -n 5

module load julia/1.9.1

MYDIR=/home/jrys/orcd/pool/metabolic-network-metrics/bacteria/results/

[ ! -d "$MYDIR/sgdMMDS/" ] && mkdir $MYDIR/sgdMMDS/

[ ! -d "$MYDIR/sgdMMDS/LGRASSMANN_ecoli_subAGORA2" ] && mkdir $MYDIR/sgdMMDS/LGRASSMANN_ecoli_subAGORA2 &&
    mkdir $MYDIR/sgdMMDS/LGRASSMANN_ecoli_subAGORA2/stress && mkdir $MYDIR/sgdMMDS/LGRASSMANN_ecoli_subAGORA2/distortion 
[ ! -d "$MYDIR/sgdMMDS/LGRASSMANN_ecoli_subAGORA2_int" ] && mkdir $MYDIR/sgdMMDS/LGRASSMANN_ecoli_subAGORA2_int &&
    mkdir $MYDIR/sgdMMDS/LGRASSMANN_ecoli_subAGORA2_int/stress && mkdir $MYDIR/sgdMMDS/LGRASSMANN_ecoli_subAGORA2_int/distortion 
[ ! -d "$MYDIR/sgdMMDS/RGRASSMANN_ecoli_subAGORA2" ] && mkdir $MYDIR/sgdMMDS/RGRASSMANN_ecoli_subAGORA2 &&
    mkdir $MYDIR/sgdMMDS/RGRASSMANN_ecoli_subAGORA2/stress && mkdir $MYDIR/sgdMMDS/RGRASSMANN_ecoli_subAGORA2/distortion
[ ! -d "$MYDIR/sgdMMDS/JACCARD_ecoli_subAGORA2" ] && mkdir $MYDIR/sgdMMDS/JACCARD_ecoli_subAGORA2 &&
    mkdir $MYDIR/sgdMMDS/JACCARD_ecoli_subAGORA2/stress && mkdir $MYDIR/sgdMMDS/JACCARD_ecoli_subAGORA2/distortion

julia xx4_Embed2.jl $SLURM_ARRAY_TASK_ID 20
