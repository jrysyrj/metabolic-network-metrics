using CSV
using DataFrames 
using Statistics
using Manifolds 
using LinearAlgebra
include("/home/jrys/sgd-MMDS/sgdMMDS.jl")

NUM_TRIALS = 5
DIMS = 1:20

savedirL    = "/home/jrys/orcd/pool/metabolic-network-metrics/bacteria/results/sgdMMDS/LGRASSMANN_ecoli_subAGORA2/"
savedirLint = "/home/jrys/orcd/pool/metabolic-network-metrics/bacteria/results/sgdMMDS/LGRASSMANN_ecoli_subAGORA2_int/"
savedirR    = "/home/jrys/orcd/pool/metabolic-network-metrics/bacteria/results/sgdMMDS/RGRASSMANN_ecoli_subAGORA2/"
savedirJ    = "/home/jrys/orcd/pool/metabolic-network-metrics/bacteria/results/sgdMMDS/JACCARD_ecoli_subAGORA2/"
loadDir     = "/home/jrys/orcd/pool/metabolic-network-metrics/bacteria/results/distances/"

N_ecoli = 389
N_subAGORA2 = 688

files = readdir(loadDir)
subAGORA_L_GRASSMANN_FILES       = [loadDir*i for i in files if !isnothing(match(r"\AsubAGORA2_LNS_\d+of\d+.csv",i))];
ecoli_L_GRASSMANN_FILES          = [loadDir*i for i in files if !isnothing(match(r"\Aecoli_LNS_\d+of\d+.csv",i))];
ecoli_subAGORA_L_GRASSMANN_FILES = [loadDir*i for i in files if !isnothing(match(r"\Aecoli_subAGORA2_LNS_\d+of\d+.csv",i))];

subAGORA_L_GRASSMANN_int_FILES       = [loadDir*i for i in files if !isnothing(match(r"\AsubAGORA2_int_LNS_\d+of\d+.csv",i))];
ecoli_L_GRASSMANN_int_FILES          = [loadDir*i for i in files if !isnothing(match(r"\Aecoli_int_LNS_\d+of\d+.csv",i))];
ecoli_subAGORA_L_GRASSMANN_int_FILES = [loadDir*i for i in files if !isnothing(match(r"\Aecoli_subAGORA2_int_LNS_\d+of\d+.csv",i))];

subAGORA_R_GRASSMANN_FILES       = [loadDir*i for i in files if !isnothing(match(r"\AsubAGORA2_RNS_\d+of\d+.csv",i))];
ecoli_R_GRASSMANN_FILES          = [loadDir*i for i in files if !isnothing(match(r"\Aecoli_RNS_\d+of\d+.csv",i))];
ecoli_subAGORA_R_GRASSMANN_FILES = [loadDir*i for i in files if !isnothing(match(r"\Aecoli_subAGORA2_RNS_\d+of\d+.csv",i))];

subAGORA_JACCARD_FILES       = [loadDir*i for i in files if !isnothing(match(r"\AsubAGORA2_jaccard_\d+of\d+.csv",i))];
ecoli_JACCARD_FILES          = [loadDir*i for i in files if !isnothing(match(r"\Aecoli_jaccard_\d+of\d+.csv",i))];
ecoli_subAGORA_JACCARD_FILES = [loadDir*i for i in files if !isnothing(match(r"\Aecoli_subAGORA2_jaccard_\d+of\d+.csv",i))];

@show length(subAGORA_JACCARD_FILES), length(ecoli_JACCARD_FILES), length(ecoli_subAGORA_JACCARD_FILES)

taskNum  = parse(Int,ARGS[1])
numTasks = parse(Int,ARGS[2])

keys = [(i,j) for i = 1:NUM_TRIALS for j ∈ DIMS]

chunks = Dict()
for i in 1:numTasks
    chunks[i] = Tuple{Int, Int}[]
end 

for i in 1:length(keys)
    push!(chunks[mod(i-1,numTasks)+1],keys[i])
end 


#########################################
###### LEFT NULLSPACE EMBEDDING 
#########################################


L_TUPLES = Tuple{Int64,Int64}[]
L_d = Float64[]

for fn in ecoli_L_GRASSMANN_FILES
    df = DataFrame(CSV.File(fn))
    global L_TUPLES = [L_TUPLES; [(df.I[k],df.J[k]) for k = 1:length(df.I)]]
    global L_d = [L_d; df.grassmann]
end
for fn in subAGORA_L_GRASSMANN_FILES
    df = DataFrame(CSV.File(fn))
    global L_TUPLES = [L_TUPLES; [(df.I[k]+N_ecoli,df.J[k]+N_ecoli) for k = 1:length(df.I)]]
    global L_d = [L_d; df.grassmann]
end
for fn in ecoli_subAGORA_L_GRASSMANN_FILES
    df = DataFrame(CSV.File(fn))
    global L_TUPLES = [L_TUPLES; [(df.I[k],df.J[k]+N_ecoli) for k = 1:length(df.I)]]
    global L_d = [L_d; df.grassmann]
end

L_weights = 1 ./ ((L_d .+ eps()) .^ 2)

for tp in chunks[taskNum]
    i = tp[1]
    k = tp[2] 

    currfnL = savedirL*lpad(k,2,"0")*"_"*lpad(i,2,"0")*".csv"

    M = Manifolds.Euclidean(k)

    if isfile(currfnL)
        continue 
    end 
   
    u, stress, distortion, exit_flag = sgdMMDS.MDS(M, L_d, L_weights, L_TUPLES, 0.1,100,100,1e-8)
    u = hcat(u...)'

    local df = DataFrame(u, :auto)
    CSV.write(currfnL, df)

    u = nothing 
    df = nothing 
    GC.gc()

    local df = DataFrame(hcat(stress), :auto)
    CSV.write(savedirL*"stress/"*lpad(k,2,"0")*"_"*lpad(i,2,"0")*".csv", df)

    stress = nothing 
    df = nothing 
    GC.gc()

    local df = DataFrame(hcat(distortion), :auto)
    CSV.write(savedirL*"distortion/"*lpad(k,2,"0")*"_"*lpad(i,2,"0")*".csv", df)

    distortion = nothing 
    df = nothing 
    GC.gc()
end 

L_d, L_weights, L_TUPLES = nothing, nothing, nothing 
GC.gc()

#########################################
###### LEFT NULLSPACE (INTERNAL) EMBEDDING 
#########################################


L_TUPLES = Tuple{Int64,Int64}[]
L_d = Float64[]

for fn in ecoli_L_GRASSMANN_int_FILES
    df = DataFrame(CSV.File(fn))
    global L_TUPLES = [L_TUPLES; [(df.I[k],df.J[k]) for k = 1:length(df.I)]]
    global L_d = [L_d; df.grassmann]
end
for fn in subAGORA_L_GRASSMANN_int_FILES
    df = DataFrame(CSV.File(fn))
    global L_TUPLES = [L_TUPLES; [(df.I[k]+N_ecoli,df.J[k]+N_ecoli) for k = 1:length(df.I)]]
    global L_d = [L_d; df.grassmann]
end
for fn in ecoli_subAGORA_L_GRASSMANN_int_FILES
    df = DataFrame(CSV.File(fn))
    global L_TUPLES = [L_TUPLES; [(df.I[k],df.J[k]+N_ecoli) for k = 1:length(df.I)]]
    global L_d = [L_d; df.grassmann]
end

L_weights = 1 ./ ((L_d .+ eps()) .^ 2)

for tp in chunks[taskNum]
    i = tp[1]
    k = tp[2] 

    currfnL = savedirLint*lpad(k,2,"0")*"_"*lpad(i,2,"0")*".csv"

    M = Manifolds.Euclidean(k)

    if isfile(currfnL)
        continue 
    end 
   
    u, stress, distortion, exit_flag = sgdMMDS.MDS(M, L_d, L_weights, L_TUPLES, 0.1,100,100,1e-8)
    u = hcat(u...)'

    local df = DataFrame(u, :auto)
    CSV.write(currfnL, df)

    u = nothing 
    df = nothing 
    GC.gc()

    local df = DataFrame(hcat(stress), :auto)
    CSV.write(savedirLint*"stress/"*lpad(k,2,"0")*"_"*lpad(i,2,"0")*".csv", df)

    stress = nothing 
    df = nothing 
    GC.gc()

    local df = DataFrame(hcat(distortion), :auto)
    CSV.write(savedirLint*"distortion/"*lpad(k,2,"0")*"_"*lpad(i,2,"0")*".csv", df)

    distortion = nothing 
    df = nothing 
    GC.gc()
end 

L_d, L_weights, L_TUPLES = nothing, nothing, nothing 
GC.gc()

#########################################
###### RIGHT NULLSPACE EMBEDDING 
#########################################

R_TUPLES = Tuple{Int64,Int64}[]
R_d = Float64[]

for fn in ecoli_R_GRASSMANN_FILES
    df = DataFrame(CSV.File(fn))
    global R_TUPLES = [R_TUPLES; [(df.I[k],df.J[k]) for k = 1:length(df.I)]]
    global R_d = [R_d; df.grassmann]
end
for fn in subAGORA_R_GRASSMANN_FILES
    df = DataFrame(CSV.File(fn))
    global R_TUPLES = [R_TUPLES; [(df.I[k]+N_ecoli,df.J[k]+N_ecoli) for k = 1:length(df.I)]]
    global R_d = [R_d; df.grassmann]
end
for fn in ecoli_subAGORA_R_GRASSMANN_FILES
    df = DataFrame(CSV.File(fn))
    global R_TUPLES = [R_TUPLES; [(df.I[k],df.J[k]+N_ecoli) for k = 1:length(df.I)]]
    global R_d = [R_d; df.grassmann]
end

R_weights = 1 ./ ((R_d .+ eps()) .^ 2)

for tp in chunks[taskNum]
    i = tp[1]
    k = tp[2] 

    currfnR = savedirR*lpad(k,2,"0")*"_"*lpad(i,2,"0")*".csv"

    M = Manifolds.Euclidean(k)

    if isfile(currfnR)
        continue 
    end 
   
    u, stress, distortion, exit_flag = sgdMMDS.MDS(M, R_d, R_weights, R_TUPLES, 0.1,30,30,1e-8)
    u = hcat(u...)'

    local df = DataFrame(u, :auto)
    CSV.write(currfnR, df)

    u = nothing 
    df = nothing 
    GC.gc()

    local df = DataFrame(hcat(stress), :auto)
    CSV.write(savedirR*"stress/"*lpad(k,2,"0")*"_"*lpad(i,2,"0")*".csv", df)

    stress = nothing 
    df = nothing 
    GC.gc()

    local df = DataFrame(hcat(distortion), :auto)
    CSV.write(savedirR*"distortion/"*lpad(k,2,"0")*"_"*lpad(i,2,"0")*".csv", df)

    distortion = nothing 
    df = nothing 
    GC.gc()
end 

R_d, R_weights, R_TUPLES = nothing, nothing, nothing 
GC.gc()


#########################################
###### JACCARD EMBEDDING 
#########################################


J_TUPLES = Tuple{Int64,Int64}[]
J_d = Float64[]

for fn in ecoli_JACCARD_FILES
    df = DataFrame(CSV.File(fn))
    global J_TUPLES = [J_TUPLES; [(df.I[k],df.J[k]) for k = 1:length(df.I)]]
    global J_d = [J_d; df.jac]
end
for fn in subAGORA_JACCARD_FILES
    df = DataFrame(CSV.File(fn))
    global J_TUPLES = [J_TUPLES; [(df.I[k]+N_ecoli,df.J[k]+N_ecoli) for k = 1:length(df.I)]]
    global J_d = [J_d; df.jac]
end
for fn in ecoli_subAGORA_JACCARD_FILES
    df = DataFrame(CSV.File(fn))
    global J_TUPLES = [J_TUPLES; [(df.I[k],df.J[k]+N_ecoli) for k = 1:length(df.I)]]
    global J_d = [J_d; df.jac]
end

J_weights = 1 ./ ((J_d .+ eps()) .^ 2)

for tp in chunks[taskNum]
    i = tp[1]
    k = tp[2] 

    currfnJ = savedirJ*lpad(k,2,"0")*"_"*lpad(i,2,"0")*".csv"

    M = Manifolds.Euclidean(k)

    if isfile(currfnJ)
        continue 
    end 
   
    u, stress, distortion, exit_flag = sgdMMDS.MDS(M, J_d, J_weights, J_TUPLES, 0.1,100,100,1e-8)
    u = hcat(u...)'

    local df = DataFrame(u, :auto)
    CSV.write(currfnJ, df)

    u = nothing 
    df = nothing 
    GC.gc()

    local df = DataFrame(hcat(stress), :auto)
    CSV.write(savedirJ*"stress/"*lpad(k,2,"0")*"_"*lpad(i,2,"0")*".csv", df)

    stress = nothing 
    df = nothing 
    GC.gc()

    local df = DataFrame(hcat(distortion), :auto)
    CSV.write(savedirJ*"distortion/"*lpad(k,2,"0")*"_"*lpad(i,2,"0")*".csv", df)

    distortion = nothing 
    df = nothing 
    GC.gc()
end 

J_d, J_weights, J_TUPLES = nothing, nothing, nothing 
GC.gc()

