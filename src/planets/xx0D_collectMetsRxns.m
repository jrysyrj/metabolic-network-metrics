clc
clear

addpath('/home/jrys/cobratoolbox')
initCobraToolbox
[solverOK,solverInstalled] = changeCobraSolver('gurobi')

basePath='/home/jrys/orcd/pool/metabolic-network-metrics/planets/results/'
fileList = dir([basePath 'consistent_planets/*_consistent.mat']);
fileNames = cell(length(fileList),1);
for k = 1:length(fileList) 
    fileNames{k} = strcat(fileList(k).folder,'/',fileList(k).name);
end 

fileNames = sort(fileNames);

allRxns = {};
allMets = {};
modelIDs = []; 


for k = 1:length(fileNames)
    fn = fileNames{k}
    load(fn);

    modelIDs = [modelIDs, {model.modelID}];

    r = sort(model.rxns);
    m = sort(model.mets);
    for j = 1:length(r) 
        if ~any(strcmp(allRxns,r{j})) 
            allRxns = [allRxns, r{j}];
        end 
    end 
    for i = 1:length(m) 
        if ~any(strcmp(allMets,m{i})) 
            allMets = [allMets, m{i}];
        end 
    end 
end 

allRxns = sort(allRxns);
allMets = sort(allMets);

info = struct();
info.rxns = allRxns;
info.mets = allMets;
info.modelIDs = modelIDs;

length(allMets) 
length(allRxns) 

save([basePath 'MetRxnInfo.mat'], '-struct', 'info');
