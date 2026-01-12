clc
clear

addpath('/home/jrys/cobratoolbox')
initCobraToolbox
[solverOK,solverInstalled] = changeCobraSolver('gurobi')

basePath='/home/jrys/orcd/pool/metabolic-network-metrics/bacteria/results/';
resultsPath=[basePath 'ecoliKOs'];
printLevel=0;

fileList = dir([resultsPath filesep 'geneGroup*.mat']);
fileNames = cell(length(fileList),1);
for k = 1:length(fileList) 
    fileNames{k} = strcat(fileList(k).folder,'/',fileList(k).name);
end 

fileNames = sort(fileNames);

isViable={};

for j = 1:length(fileNames)
    fn = fileNames{j}
    load(fn);
    sol = optimizeCbModel(model);
    isViable{j}(1) = sol.f > 0;
end

save([resultsPath filesep 'Model_viability'],'isViable','-v7.3');
