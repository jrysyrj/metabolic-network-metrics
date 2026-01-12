using LinearAlgebra, CSV, DataFrames

N_ecoli = 389

baseDir = "/home/jrys/orcd/pool/metabolic-network-metrics/bacteria/results/"
loadDir_ecoli_L = baseDir*"ecoli_LNS/"
loadDir_ecoli_int_L = baseDir*"ecoli_int_LNS/"
loadDir_ecoli_R = baseDir*"ecoli_RNS/"

L_nullity = zeros(Float64,N_ecoli)
L_int_nullity = zeros(Float64,N_ecoli)
R_nullity = zeros(Float64,N_ecoli)

for i = 1:N_ecoli 
    @show i
    fn = loadDir_ecoli_L*lpad(i,3,"0")*".csv"
    df = DataFrames.DataFrame(CSV.File(fn,header=0))
    S = Matrix(df)
    L_nullity[i] = size(S,2)
    
    fn, df, S = nothing, nothing, nothing  
    GC.gc()
    
    fn = loadDir_ecoli_int_L*lpad(i,3,"0")*".csv"
    df = DataFrames.DataFrame(CSV.File(fn,header=0))
    S = Matrix(df)
    L_int_nullity[i] = size(S,2)
    
    fn, df, S = nothing, nothing, nothing  
    GC.gc()

    fn = loadDir_ecoli_R*lpad(i,3,"0")*".csv"
    df = DataFrames.DataFrame(CSV.File(fn,header=0))
    S = Matrix(df)
    R_nullity[i] = size(S,2)
    
    fn, df, S = nothing, nothing, nothing  
    GC.gc()
end 

df = DataFrames.DataFrame(Lnullity=L_nullity, Lintnullity=L_int_nullity, Rnullity=R_nullity)
CSV.write(baseDir*"ecoli_nullity.csv",df)


N_subAGORA2 = 688

baseDir = "/home/jrys/orcd/pool/metabolic-network-metrics/bacteria/results/"
loadDir_subAGORA2_L = baseDir*"subAGORA2_LNS/"
loadDir_subAGORA2_int_L = baseDir*"subAGORA2_int_LNS/"
loadDir_subAGORA2_R = baseDir*"subAGORA2_RNS/"

L_nullity = zeros(Float64,N_subAGORA2)
L_int_nullity = zeros(Float64,N_subAGORA2)
R_nullity = zeros(Float64,N_subAGORA2)

for i = 1:N_subAGORA2 
    @show i
    fn = loadDir_subAGORA2_L*lpad(i,3,"0")*".csv"
    df = DataFrames.DataFrame(CSV.File(fn,header=0))
    S = Matrix(df)
    L_nullity[i] = size(S,2)
    
    fn, df, S = nothing, nothing, nothing  
    GC.gc()
    
    fn = loadDir_subAGORA2_int_L*lpad(i,3,"0")*".csv"
    df = DataFrames.DataFrame(CSV.File(fn,header=0))
    S = Matrix(df)
    L_int_nullity[i] = size(S,2)
    
    fn, df, S = nothing, nothing, nothing  
    GC.gc()

    fn = loadDir_subAGORA2_R*lpad(i,3,"0")*".csv"
    df = DataFrames.DataFrame(CSV.File(fn,header=0))
    S = Matrix(df)
    R_nullity[i] = size(S,2)
    
    fn, df, S = nothing, nothing, nothing  
    GC.gc()
end 

df = DataFrames.DataFrame(Lnullity=L_nullity, Lintnullity=L_int_nullity, Rnullity=R_nullity)
CSV.write(baseDir*"subAGORA2_nullity.csv",df)
