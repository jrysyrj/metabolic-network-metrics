using LinearAlgebra, CSV, DataFrames, Glob


taskNum  = parse(Int,ARGS[1])
numTasks = parse(Int,ARGS[2])

loadDir = nothing 
try 
    global loadDir = ARGS[3]
catch
    global loadDir  = "~/"
end
loadDir2 = nothing 
try
    global loadDir2 = ARGS[4] 
catch
    global loadDir2 = "~/"
end

savePrefix = nothing
try
    global savePrefix = ARGS[5]
catch
    global savePrefix = "~/"
end

println(loadDir)
println(loadDir2)
println(savePrefix)

files = sort(glob("*.csv",loadDir2))
N = length(files)

keys = [(j,i) for i in 1:N for j in 1:i-1]

chunks = Dict()
for i in 1:numTasks
    chunks[i] = Tuple{Int, Int}[]
end 
for i in 1:length(keys)
    push!(chunks[mod(i-1,numTasks)+1],keys[i])
end 

keys = chunks[taskNum]

I = [keys[i][1] for i=1:length(keys)]
J = [keys[i][2] for i=1:length(keys)]

g = zeros(Float64,length(keys))
κ = zeros(Float64,length(keys))
ρ = zeros(Float64,length(keys))
πdist = zeros(Float64,length(keys))
ϵ = zeros(Float64,length(keys))

for i = 1:length(keys)
    k = keys[i]
    fn = loadDir*lpad(k[1],3,"0")*"-"*lpad(k[2],3,"0")*".csv"
    df = DataFrames.DataFrame(CSV.File(fn,header=0))

    fnS = loadDir2*lpad(k[1],3,"0")*".csv"
    fnT = loadDir2*lpad(k[2],3,"0")*".csv"
    dfS = DataFrame(CSV.File(fnS, header=0))
    dfT = DataFrame(CSV.File(fnT, header=0))
    S = Matrix(dfS)
    T = Matrix(dfT)

    abs_diff_dim = abs(size(S,2) - size(T,2))
    fnS, fnT, dfS, dfT = nothing, nothing, nothing, nothing  
    GC.gc()

    θ = df[!,1]
    sinθ = clamp.(sin.(θ),0.0,1.0)
    sinθdiv2 = clamp.(sin.(θ/2),0.0,1.0)

    ϵ[i] = sqrt(abs_diff_dim)
    g[i] = sqrt(sum(θ .^ 2) + abs_diff_dim * (π^2 / 4))
    κ[i] = sqrt(sum(sinθ .^ 2) + abs_diff_dim)
    ρ[i] = sqrt(2*sum(sinθdiv2 .^ 2) + abs_diff_dim)
    if iszero(abs_diff_dim)
        πdist[i] = clamp(sin(maximum(θ)),0.0,1.0)
    else 
        πdist[i] = 1
    end
    fn, df, θ  = nothing, nothing, nothing 
    cosθ, sinθ, sinθdiv2 = nothing, nothing, nothing 
    GC.gc()
end 

df = DataFrames.DataFrame(I=I, J=J, grassmann=g, chordal=κ, procrustes=ρ, projection=πdist, dimgap=ϵ)

CSV.write(savePrefix*"_"*string(taskNum)*"of"*string(numTasks)*".csv",df)
