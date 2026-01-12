using CSV, DataFrames 

fn = "/home/jrys/orcd/pool/metabolic-network-metrics/bacteria/results/subAGORA2_rxns_present.csv"
df = DataFrame(CSV.File(fn,header=0))
M = Matrix(df)
N = size(M,1)

keys = [(j,i) for i in 1:N for j in 1:i-1]
I = [keys[i][1] for i=1:length(keys)]
J = [keys[i][2] for i=1:length(keys)]

jac = zeros(Float64,length(keys))

for i = 1:length(keys)
    @show k = keys[i]
    b = M[k[1],:]
    c = M[k[2],:] 
    jac[i] = 1 - (sum(b .& c)/sum(b .| c))
end 

df = DataFrames.DataFrame(I=I, J=J, jac=jac)

CSV.write("/home/jrys/orcd/pool/metabolic-network-metrics/bacteria/results/distances/subAGORA2_jaccard_1of1.csv",df)

#######

fn = "/home/jrys/orcd/pool/metabolic-network-metrics/bacteria/results/subAGORA2_rxns_present_int.csv"
df = DataFrame(CSV.File(fn,header=0))
M = Matrix(df)
N = size(M,1)

keys = [(j,i) for i in 1:N for j in 1:i-1]
I = [keys[i][1] for i=1:length(keys)]
J = [keys[i][2] for i=1:length(keys)]

jac = zeros(Float64,length(keys))

for i = 1:length(keys)
    @show k = keys[i]
    b = M[k[1],:]
    c = M[k[2],:] 
    jac[i] = 1 - (sum(b .& c)/sum(b .| c))
end 

df = DataFrames.DataFrame(I=I, J=J, jac=jac)

CSV.write("/home/jrys/orcd/pool/metabolic-network-metrics/bacteria/results/distances/subAGORA2_int_jaccard_1of1.csv",df)