using CSV, DataFrames, LinearAlgebra, SparseArrays, MAT

### --- BANK --- ###

planets = ["ModernEarthY&D99", "Titan2016", "ModernMars", "ModernEarthSanders", 
            "Venus", "EarlyMars", "Jupiter", "Pluto"]
remove_species = ["RAIN","VSTRDUST","MARSDUST","RAYEARTH","RAYCO2","VENRAY","HV","PROD","U","M","ADDUCT",
                    "SGA","JDUST","THERM","TDUST1","TDUST2","TDUST3","TDUST4","TDUST5","TDUST6","TDUST7",
                    "TDUST8","TDUST8b","RAYJUP","V","VHAZE1","VHAZE2","W","X","Y"]
    

function parse_species(str)
    ans = ""
    regex = r"(\d+)\)\s*([\w\d\(\)\+\-\,]*)\s*\d*.*"
    # Use a regular expression to extract the number, molecule, and state
    match = eachmatch(regex, str) |> first
    if match !== nothing
        number = parse(Int, match.captures[1])  # Extract and convert the number
        molecule = match.captures[2]           # Extract the molecule
        ans = molecule
    else
        error("Failed to parse species: $str")
    end
    ans = replace(ans,"CL"=>"Cl", count=1)
    ans = replace(ans,"BR"=>"Br", count=1)
    return ans
end

function parse_reaction(str)
    if str[1:3] == "koo"
        return 
    end
    regex = r"(\d+\))?\s*([\w\s\d\(\)\+\-\,]*)\s*=\s*([\w\s\d\(\)\+\-\,]*)\s*((?:k.*|\(.*))?$"
    stoich_regex = r"(\d*\s*)([A-Za-z0-9\+\-]+)"
    ans = (reactants=[],products=[])
    # Match the reaction pattern
    match = eachmatch(regex, str) |> first
    if match !== nothing
        # Extract matched reactants and products strings
        reactants_str = match.captures[2]
        products_str = match.captures[3]

        # Split reactants and products
        reactants = string.([strip(r) for r in split(reactants_str, r"\s\+\s")]) # Split by "+" and remove spaces
        products = string.([strip(r) for r in split(products_str, r"\s\+\s")])    # Split by "+" and remove spaces

        reactants_p = [] 
        for m in reactants 
            stoich = 1
            m̃ = ""
            foo = split(m,r"\s")
            if length(foo) == 1
                m̃ = foo[1]
            else
                stoich = parse(Int,foo[1])
                m̃ = foo[2]
            end
            m̃ = replace(m̃,"CL"=>"Cl", count=1)
            m̃ = replace(m̃,"BR"=>"Br", count=1)
            push!(reactants_p,(stoich=stoich, species=m̃))
        end

        products_p = [] 
        for m in products 
            stoich = 1
            m̃ = ""
            foo = split(m,r"\s")
            if length(foo) == 1
                m̃ = foo[1]
            else
                stoich = parse(Int,foo[1])
                m̃ = foo[2]
            end
            m̃ = replace(m̃,"CL"=>"Cl", count=1)
            m̃ = replace(m̃,"BR"=>"Br", count=1)
            push!(products_p,(stoich=stoich, species=m̃))
        end

        filter!(e->e.species∉remove_species,reactants_p) 
        filter!(e->e.species∉remove_species,products_p) 
        
        reactants = reactants_p
        products = products_p

        if length(reactants) == 1 && reactants[1] == ""
            println("Missing reactants information: $str") 
            return 
        end
        if length(reactants) > 1 && length(products) == 0
            println("Cannot combine species and produce nothing!: $str") 
            return 
        end
        if length(products) > 1 && length(reactants) == 0
            println("Cannot decompose species from nothing!: $str") 
            return 
        end
        if length(products) == 1 && products[1] == ""
            println("Missing products information: $str") 
            return 
        end
        if reactants == products
            println("Assuming reactants = products corresponds to exchange: $str")
            products = typeof(reactants[1])[]
        end

        reactant_species = [e.species for e in reactants]
        product_species = [e.species for e in products]

        shared = string.(collect(intersect(Set(reactant_species),Set(product_species))))
        if length(shared) != 0
            for k in shared
                idx_r = findall(reactant_species .== k) |> first
                idx_p = findall(product_species .== k)  |> first
                ΔSTOICH = reactants[idx_r].stoich -  products[idx_p].stoich  
                if ΔSTOICH  == 0 
                    deleteat!(reactants,idx_r)
                    deleteat!(products,idx_p)
                elseif ΔSTOICH  > 0
                    foo = (stoich = ΔSTOICH, species = reactants[idx_r].species)
                    reactants[idx_r] = foo
                    deleteat!(products,idx_p)
                elseif ΔSTOICH  < 0
                    foo = (stoich = -ΔSTOICH, species = products[idx_p].species)
                    products[idx_p] = foo
                    deleteat!(reactants,idx_r)
                end 
            end 
        end 
        
        
        ans = (reactants=reactants, products=products)
    else
        error("Failed to parse reaction: $str")
    end
    return ans
end

loadDir = "/home/jrys/orcd/pool/metabolic-network-metrics/planets/data/"
saveDir = "/home/jrys/orcd/pool/metabolic-network-metrics/planets/results/networks_MAT/"

for planet_name in planets 
    fn = loadDir*planet_name*".edited.out"

    species = []
    reactions = []

    open(fn) do f
        @show f 
        while ! eof(f)    
            s = strip(readline(f))
            if s == "SPECIES:"
                while true 
                    try 
                        push!(species,parse_species(strip(readline(f))))
                    catch
                        break 
                    end
                end
            end
            if s == "REACTIONS:"
                while true 
                    try 
                        ans = parse_reaction(strip(readline(f)))
                        if ans != nothing 
                            push!(reactions,ans)
                        end
                    catch
                        break 
                    end
                end

            end
        end
    end

    filter!(e->e∉remove_species,species) 

    @show length(species) 
    @show length(reactions) 

    mets = sort(species)
    rxns = [join(k.reactants," + ")*" = "*join(k.products," + ") for k in reactions]
    sp = sortperm(rxns)
    rxns = rxns[sp]
    reactions = reactions[sp] 

    stoich_regex = r"(\d*\s*)([A-Za-z0-9\+\-]+)"

    STOICH_MATRIX = spzeros(length(mets),length(rxns))
    for (idr, k) in enumerate(reactions)
        participants = []
        for m in k.reactants
            push!(participants,(-m.stoich,m.species))
        end
        for m in k.products
            push!(participants,(m.stoich,m.species))
        end
        for (stoich,sp) in participants
            idm = findall(mets .== sp)[1]
            STOICH_MATRIX[idm,idr] = stoich 
        end 
    end

    model = Dict("mets"=>mets,"S"=>STOICH_MATRIX,"modelID"=>planet_name)
    matwrite(saveDir*planet_name*".mat",Dict("model"=>model))
end