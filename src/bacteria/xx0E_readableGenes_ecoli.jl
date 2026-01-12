using MAT
using CSV, DataFrames

function correct_PEG_ID(s::AbstractString)
    gene_parts = split(s, ".")
    gene_parts = [replace(p, "\"" => "") for p in gene_parts if p != "peg"]
    last_part = gene_parts[end]
    underscore_idx = findfirst(==('_'), last_part)
    if underscore_idx !== nothing
        last_part = last_part[1:underscore_idx-1]
    end
    gene_id = string(gene_parts[1], ".", gene_parts[2], ".peg.", last_part)
    return gene_id
end


vars = matread("/home/jrys/orcd/pool/metabolic-network-metrics/bacteria/results/Escherichia_coli_str_K_12_substr_MG1655_consistent.mat")

AGORA_GENE_NAMES = vcat(vars["model"]["genes"]...)

PEG_NAMES = String[]
nonPEG_NAMES = String[]

for gene in AGORA_GENE_NAMES
    if occursin("peg",gene) push!(PEG_NAMES,gene)
    else push!(nonPEG_NAMES,gene)
    end
end

println(length(PEG_NAMES),", ",length(nonPEG_NAMES))

correctPEG_NAMES = String[]
for gene in PEG_NAMES
    push!(correctPEG_NAMES,correct_PEG_ID(gene))
end 


df_BVBRC = DataFrame(CSV.File("/home/jrys/orcd/pool/metabolic-network-metrics/bacteria/data/ecoli_BVBRC_genome_feature.csv"))
BRC_ID = [string(k[5:end]) for k in df_BVBRC[!,4]]
GENE_ID = df_BVBRC[!,19]

FOUND_IN_BRC = String[]
BRC_MAPPED = Union{String,Missing}[]
notFOUND_IN_BRC = String[]

for gene in correctPEG_NAMES
    if gene in BRC_ID
        idx = findfirst(==(gene),BRC_ID)
        push!(FOUND_IN_BRC,gene)
        push!(BRC_MAPPED,GENE_ID[idx])
    else push!(notFOUND_IN_BRC,gene)
    end
end

println(length(FOUND_IN_BRC),", ",length(notFOUND_IN_BRC))


SIMPLE_TRANSLATED = fill("",length(AGORA_GENE_NAMES))

for k = 1:length(AGORA_GENE_NAMES)
    gene = AGORA_GENE_NAMES[k]
    result = gene
    if occursin("peg",gene) 
        gene = correct_PEG_ID(gene)
        result = gene
        idx = findfirst(==(gene),BRC_ID)
        if !ismissing(GENE_ID[idx])
            result = GENE_ID[idx]
        end
    end
    SIMPLE_TRANSLATED[k] = result 
    println(AGORA_GENE_NAMES[k], "\t", result)
end

df = DataFrames.DataFrame(AGORA=AGORA_GENE_NAMES, READABLE=SIMPLE_TRANSLATED)
CSV.write("/home/jrys/orcd/pool/metabolic-network-metrics/bacteria/results/SIMPLE_MATCH_ecoli_genes.csv",df)


