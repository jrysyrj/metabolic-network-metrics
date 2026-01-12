#!/bin/bash

#SBATCH -o jaccard.sh.log-%j

module load julia/1.9.1

julia xx0K_jaccard_ecoli.jl 
julia xx0K_jaccard_subAGORA2.jl 
julia xx0K_jaccard_ecoli_subAGORA2.jl 