clc
clear

addpath('/home/jrys/cobratoolbox')
initCobraToolbox
[solverOK,solverInstalled] = changeCobraSolver('gurobi')

basePath='/home/jrys/orcd/pool/metabolic-network-metrics/bacteria'
resultsPath=[basePath '/results/consistent_subAGORA2']
mkdir(resultsPath)
printLevel=0;

fileList = dir([basePath '/data/subAGORA2/*.mat'])
fileNames = cell(length(fileList),1);
for k = 1:length(fileList) 
    fileNames{k} = strcat(fileList(k).folder,'/',fileList(k).name);
end 
fileNames = sort(fileNames);

% param for consistency checks 
param=struct;
feasTol = getCobraSolverParams('LP', 'feasTol');
param.('feasTol')= feasTol;
param.('epsilon')=feasTol*100;
param.('modeFlag')=0;
param.('method')='fastcc';

% flags for rxn removal 
irrevFlag=0;
metFlag=1;

massImbalancedRxns = [];
massImbalancedRxnFormulas = [];

for j = 1:length(fileNames)
    fn = fileNames{j}
    load(fn);
    
    [model, removedRxnInd, keptRxnInd] = checkDuplicateRxn(model);
    
    
    [massImbalance, imBalancedMass, imBalancedCharge, imBalancedRxnBool, elements, missingFormulaeBool, balancedMetBool]  = checkMassChargeBalance(model);
    model = findSExRxnInd(model,[],printLevel-1);
    exRxns = find(~model.SIntRxnBool);
    
    massImbalancep = massImbalance;
    massImbalancep(exRxns,:) = [];
    imBalancedMassp = imBalancedMass;
    imBalancedMassp(exRxns) = [];
    modelrxnsp = model.rxns;
    modelrxnsp(exRxns) = [];
    idxs = find(~strcmp(imBalancedMassp,''));
    idxs = idxs(~contains(imBalancedMassp(idxs),'X'));
    
    exclude = {'TECAAE','TECAGE','TECAUE'};
    keep = [];
    for i = 1:length(idxs)
        if ~contains(exclude,modelrxnsp(idxs(i)))
            keep = [keep; i];
        end
    end
    idxs = idxs(keep);
    potentially_mass_imbalanced = modelrxnsp(idxs)
    massImbalancedRxns = [massImbalancedRxns; potentially_mass_imbalanced];
    massImbalancedRxnFormulas = [massImbalancedRxnFormulas; printRxnFormula(model,potentially_mass_imbalanced)];
end 

massImbalancedRxns_u = unique(massImbalancedRxns,'stable');
massImbalancedRxnFormulas_u = unique(massImbalancedRxnFormulas,'stable');

table={'Reaction name','Reaction formula'};
table(2:length(massImbalancedRxns_u)+1,1)=massImbalancedRxns_u;
table(2:length(massImbalancedRxns_u)+1,2)=massImbalancedRxnFormulas_u;

cell2csv(['MassImbalanced.csv'],table)