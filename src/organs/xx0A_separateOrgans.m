clc
clear

%addpath('/home/jrys/cobratoolbox')
%initCobraToolbox
%[solverOK,solverInstalled] = changeCobraSolver('gurobi')

basePath='/home/jrys/orcd/pool/metabolic-network-metrics/organs'
resultsPath=[basePath '/results/organ_recons']
mkdir(resultsPath)
printLevel=0;

fileList = dir([basePath '/data/OrganAtlas_*.mat'])
fileNames = cell(length(fileList),1);
for k = 1:length(fileList)
    fileNames{k} = strcat(fileList(k).folder,'/',fileList(k).name);
end
fileNames = sort(fileNames)

for k = 1:length(fileNames)
    fn = fileNames{k}
    data = load(fn);
    fieldNames = fieldnames(data);
    idx = find(startsWith(fieldNames,'OrganCompendium'),1);
    desiredFieldName = fieldNames{idx};
    OrganCompendium = data.(desiredFieldName);
    which_organs = fieldnames(OrganCompendium);
    which_organs = which_organs(1:end-2);
    sex = OrganCompendium.sex;
    for j = 1:length(which_organs)
        which_organs{j};
        model = OrganCompendium.(which_organs{j}).modelAllComp;	    
        model.modelID = [sex '_' which_organs{j}];
        resultsFileName=[resultsPath filesep model.modelID '.mat']
        save(resultsFileName,'model');
    end
end

