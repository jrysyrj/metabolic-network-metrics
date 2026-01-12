clc
clear

addpath('/home/jrys/cobratoolbox')
initCobraToolbox
[solverOK,solverInstalled] = changeCobraSolver('gurobi')

basePath='/home/jrys/orcd/pool/metabolic-network-metrics/bacteria/results/'

model = readCbModel([basePath 'Escherichia_coli_str_K_12_substr_MG1655_consistent.mat'])
df = readtable([basePath 'SIMPLE_MATCH_ecoli_genes.csv'], 'TextType', 'string');
map = containers.Map(df{:,1}, df{:,2});

MODEL_GENES = model.genes;

absent_GENES = string.empty;
for i = 1:length(MODEL_GENES)
    gene = string(MODEL_GENES{i});
    if ~isKey(map, gene)
        absent_GENES(end+1) = gene;
    end
end
disp("Number of absent genes: " + length(absent_GENES));

for k = 1:length(MODEL_GENES)
    gene = string(MODEL_GENES{k});  % convert to string if needed
    if ~any(absent_GENES == gene)
        MODEL_GENES{k} = char(map(gene));
    end
end

MODEL_GRRULES = model.grRules;  % Cell array of strings
rxnGeneMat = model.rxnGeneMat;  % Sparse or logical matrix

for k = 1:length(MODEL_GRRULES)
    foo = string(MODEL_GRRULES{k});  % Convert to string for easier replacement

    % Find indices of genes involved in reaction k
    idx = find(rxnGeneMat(k,:));  % Returns a row vector

    if isempty(idx)
        continue;
    end

    for i = idx
        gene = string(model.genes{i});  % Gene name as string

        if ~any(absent_GENES == gene)
            foo = strrep(foo, gene, map(gene));
        end
    end
    MODEL_GRRULES{k} = char(foo);
end

model.genes = MODEL_GENES;
model.grRules = MODEL_GRRULES;

model.genes = unique(model.genes);
model = rmfield(model, 'rules');
model = rmfield(model, 'rxnGeneMat');
model = generateRules(model);
model = buildRxnGeneMat(model);

writeCbModel(model,'format','mat','fileName',[basePath 'ECOLI_TRANSLATED.mat'])
