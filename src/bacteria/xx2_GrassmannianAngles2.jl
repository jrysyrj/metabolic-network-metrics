using LinearAlgebra, CSV, DataFrames, Glob

function principal_angles(fnS,fnT)
    dfS = DataFrame(CSV.File(fnS, header=0))
    dfT = DataFrame(CSV.File(fnT, header=0))

    S = Matrix(dfS)
    T = Matrix(dfT)

    fnS, fnT, dfS, dfT = nothing, nothing, nothing, nothing  
    GC.gc()

    foo = S' * T
    σ = svdvals(foo)
    σ = clamp.(σ,0.0,1.0)
    σ[isapprox.(σ, 1.0)] .= 1.0 
    σ[isapprox.(σ, 0.0)] .= 0.0
    θ = acos.(σ)

    S, T = nothing, nothing, nothing, nothing  
    foo = nothing 
    GC.gc()
    return θ
end 


taskNum  = parse(Int,ARGS[1])
numTasks = parse(Int,ARGS[2])

loadDir1 = nothing 
try 
    global loadDir1 = ARGS[3] 
catch
    global loadDir1 = "~/"
end
println(loadDir1)

loadDir2 = nothing 
try 
    global loadDir2 = ARGS[4] 
catch
    global loadDir2 = "~/"
end
println(loadDir2)

saveDir = nothing 
try 
    global saveDir = ARGS[5] 
catch
    global saveDir = "~/"
end
println(saveDir)


files1 = sort(glob("*.csv",loadDir1))
N1 = length(files1)
files2 = sort(glob("*.csv",loadDir2))
N2 = length(files2)
keys = [(i,j) for i in 1:N1 for j in 1:N2]

chunks = Dict()
for i in 1:numTasks
    chunks[i] = Tuple{Int, Int}[]
end 

for i in 1:length(keys)
    push!(chunks[mod(i-1,numTasks)+1],keys[i])
end 

for i in chunks[taskNum]
    fn = saveDir*lpad(i[1],3,"0")*"-"*lpad(i[2],3,"0")*".csv"
    if isfile(fn)
        continue 
    end 
    df = DataFrames.DataFrame(angles = principal_angles(files1[i[1]],files2[i[2]]))
    CSV.write(fn,df,header=false)
end

