using CSV
using DataFrames 
using Statistics
using Manifolds 
using LinearAlgebra
using PhyloNetworks
include("/home/jrys/sgd-MMDS/sgdMMDS.jl")

NUM_TRIALS = 5
DIMS = 1:20

savedir = "/home/jrys/orcd/pool/metabolic-network-metrics/bacteria/results/sgdMMDS/PHYLOGENETIC_subAGORA2/"

tree_file = "/home/jrys/orcd/pool/metabolic-network-metrics/bacteria/data/subAGORA2_genomic_outgroup_rooted.tre.treefile"
T = PhyloNetworks.readMultiTopology(tree_file)[1]
maxlikelihood = sqrt.(pairwisetaxondistancematrix(T; keepInternal=false))

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

tri = triu!(trues(size(maxlikelihood)),1)
CartInds = findall(tri)

TUPLES = [(k[1],k[2]) for k in CartInds]
d = maxlikelihood[CartInds]

weights = 1 ./ ((d .+ eps()) .^ 2)

for tp in chunks[taskNum]
    i = tp[1]
    k = tp[2] 

    currfn = savedir*lpad(k,2,"0")*"_"*lpad(i,2,"0")*".csv"

    M = Manifolds.Euclidean(k)

    if isfile(currfn)
        continue 
    end 
   
    u, stress, distortion, exit_flag = sgdMMDS.MDS(M, d, weights, TUPLES, 0.1,100,100,1e-8)
    u = hcat(u...)'

    local df = DataFrame(u, :auto)
    CSV.write(currfn, df)

    u = nothing 
    df = nothing 
    GC.gc()

    local df = DataFrame(hcat(stress), :auto)
    CSV.write(savedir*"stress/"*lpad(k,2,"0")*"_"*lpad(i,2,"0")*".csv", df)

    stress = nothing 
    df = nothing 
    GC.gc()

    local df = DataFrame(hcat(distortion), :auto)
    CSV.write(savedir*"distortion/"*lpad(k,2,"0")*"_"*lpad(i,2,"0")*".csv", df)

    distortion = nothing 
    df = nothing 
    GC.gc()
end 

