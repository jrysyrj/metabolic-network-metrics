using CSV, DataFrames, SparseArrays, LinearAlgebra, MAT, Glob

function my_nullspace(A::AbstractVecOrMat; atol::Real = 0.0, rtol::Real = (min(size(A, 1), size(A, 2))*eps(real(float(oneunit(eltype(A))))))*iszero(atol), perturb::Int=0)
    m, n = size(A, 1), size(A, 2)
    (m == 0 || n == 0) && return Matrix{eigtype(eltype(A))}(I, n, n)
    SVD = svd(A; full=true,alg=LinearAlgebra.QRIteration())
    @show tol = max(atol, SVD.S[1]*rtol)
    indstart = clamp(sum(s -> s .> tol, SVD.S) + 1 + perturb, 1, size(SVD.Vt,1))
    return copy((@view SVD.Vt[indstart:end,:])')
end

function KernelFromFile(fnS::String,isRight::Bool, nmets::Int, nrxns::Int, Δ::Int)
    dfS = DataFrame(CSV.File(fnS, header=0))
    I = dfS.Column1
    J = dfS.Column2
    val = dfS.Column3
    S = Array(sparse(I,J,val,nmets,nrxns))
    
    # DO THIS FOR LEFT NULLSPACE
    if !isRight
        S = S'
    end

    WHERE_EMPTY = [all(iszero.(S[:,i])) for i = 1:size(S,2)]
    S̃ = S[:,.!WHERE_EMPTY]
    WHERE_EMPTY_2 = [all(iszero.(S̃[i,:])) for i = 1:size(S̃,1)]
    S̃ = S̃[.!WHERE_EMPTY_2,:]

    try 
        A = my_nullspace(S̃;perturb=Δ)
        Ã = zeros(size(S,2),size(A,2))
        Ã[.!WHERE_EMPTY,:] .= A 
        return Ã
    catch e 
        println("Problem with "*fnS)
    end 
    return Float64[]
end 

taskNum  = parse(Int,ARGS[1])
numTasks = parse(Int,ARGS[2])

isRight = nothing 
try 
    global isRight = parse(Bool,ARGS[3]) 
catch
    global isRight = true 
end

loadDir = nothing
try
    global loadDir = ARGS[4]
catch
    global loadDir = "~/"
end
println(loadDir)

saveDir = nothing 
try 
    global saveDir = ARGS[5] 
catch
    global saveDir = "~/"
end
println(saveDir)

nmets = nothing 
try 
    global nmets = parse(Int,ARGS[6])
catch
    global nmets = 0
end
println(nmets)

nrxns = nothing 
try 
    global nrxns = parse(Int,ARGS[7])
catch
    global nrxns = 0
end
println(nrxns)

Δs = [-2, -1, 1, 2]
chunks = Dict()

for i in 1:numTasks
    chunks[i] = Tuple[]
end 

files = sort(glob("*.txt",loadDir))
N = length(files)
for j = 1:length(Δs)
    for i = 1:N 
        push!(chunks[mod(i-1,numTasks)+1],(i, j))
    end 
end

for (i, j) in chunks[taskNum] 
    fn = saveDir*string(j)*"/"*lpad(i,3,"0")*".csv" 

    A = KernelFromFile(files[i], isRight, nmets, nrxns, Δs[j])
    if !isempty(A) 
        if isRight 
            @show size(A), nrxns
            @assert size(A,1) == nrxns
        else
            @show size(A), nmets
            @assert size(A,1) == nmets
        end
    end
    df = DataFrames.DataFrame(A, :auto)
    CSV.write(fn,df,header=false)
    A, df = nothing, nothing
    GC.gc()
end 