#!/bin/bash

#SBATCH -o embed.sh.log-%j-%a
#SBATCH -a 1-19
#SBATCH -n 4

module load julia/1.9.1

MYDIR=/home/jrys/orcd/pool/metabolic-network-metrics/bacteria/results/

[ ! -d "$MYDIR/sgdMMDS/" ] && mkdir $MYDIR/sgdMMDS/
[ ! -d "$MYDIR/sgdMMDS/PHYLOGENETIC_subAGORA2" ] && mkdir $MYDIR/sgdMMDS/PHYLOGENETIC_subAGORA2
[ ! -d "$MYDIR/sgdMMDS/PHYLOGENETIC_subAGORA2/stress" ] && mkdir $MYDIR/sgdMMDS/PHYLOGENETIC_subAGORA2/stress
[ ! -d "$MYDIR/sgdMMDS/PHYLOGENETIC_subAGORA2/distortion" ] && mkdir $MYDIR/sgdMMDS/PHYLOGENETIC_subAGORA2/distortion

julia xx4_Embed_Phylogenetic.jl $SLURM_ARRAY_TASK_ID 19

