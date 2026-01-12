clc
clear

addpath('/home/jrys/cobratoolbox')
initCobraToolbox
[solverOK,solverInstalled] = changeCobraSolver('gurobi')

basePath='/home/jrys/orcd/pool/metabolic-network-metrics/bacteria'
resultsPath=[basePath '/results/']
mkdir(resultsPath)
printLevel=0;

fileList = dir([basePath '/data/*.mat'])
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
    
    if any(strcmp(model.rxns,'N2Ormqpp'))
        model = addReaction(model, 'N2Ormqpp',...
            'metaboliteList', {'mql8[c]', 'n2o[c]', 'h2o[c]', 'mqn8[c]','n2[c]'},...
            'stoichCoeffList', [-1; -1; 1; 1; 1]);
    end
    if any(strcmp(model.rxns,'MPPP9MMEOR3'))
        model = addReaction(model, 'MPPP9MMEOR3',...
            'metaboliteList', {'nadph[c]', 'o2[c]', 'oxomppp9mme[c]', 'dvprotochloroph[c]','h2o[c]', 'nadp[c]'},...
            'stoichCoeffList', [-1; -1; -1; 1; 2; 1]);
    end

    if any(strcmp(model.rxns,'GLCS3'))
        model = addReaction(model, 'GLCS3',...
            'metaboliteList', {'adpglc[c]', 'glycogenb[c]', 'adp[c]', 'glycogen[c]','h[c]'},...
            'stoichCoeffList', [-1; -1; 1; 1; 1]);
    end

    if any(strcmp(model.rxns,'GLCP3'))
        model = addReaction(model, 'GLCP3',...
            'metaboliteList', {'pi[c]', 'glycogen[c]', 'g1p[c]', 'glycogenb[c]'},...
            'stoichCoeffList', [-1; -1; 1; 1]);
    end
    
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
    
    exclude = {'TECAAE','TECAGE','TECAUE','LPCHOL23AM','HMR_9498','GLYCPCHOLAH','MAN1PT2','PMANMi'};
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

cell2csv(['MassImbalanced_ecoli_bsub.csv'],table)
