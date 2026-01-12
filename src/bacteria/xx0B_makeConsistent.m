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

    
tmpStats={};
tmpATP={};
tmpSFConsist={};


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
    mass_imbalanced = modelrxnsp(idxs);
    model = removeRxns(model,mass_imbalanced,irrevFlag,metFlag);

    tmpStats{j}(1)=length(model.rxns);
    tmpStats{j}(2)=length(model.mets);
    tmpStats{j}(3)=length(model.genes);
    
    [atpFluxAerobic, atpFluxAnaerobic] = testATP(model); % mmol/gDW/h
    tmpATP{j}(1)=atpFluxAerobic;
    tmpATP{j}(2)=atpFluxAnaerobic;
    
    modelOri = model;
    model = findSExRxnInd(model,[],printLevel-1);
    exRxns = find(~model.SIntRxnBool);
    
    % remove constraints on exchange reactions
    model.lb(exRxns,1)=-1000;
    model.ub(exRxns,1)=1000;
    
    % remove forced flux constraints
    model.lb(model.lb>0)=0;

    % remove bounds on ATP production
    model=changeRxnBounds(model,'DM_atp_c_',1000,'u');

    % test stochiometric and flux consistency
    [fluxConsistentMetBool, fluxConsistentRxnBool, fluxInConsistentMetBool, fluxInConsistentRxnBool] = findFluxConsistentSubset(model,param);
    % exclude exchange and demand reactions
    model = findSExRxnInd(model,[],printLevel-1);
    exRxns = find(~model.SIntRxnBool);
    
    foo = fluxConsistentRxnBool;
    bar = fluxInConsistentRxnBool;
    foo(exRxns) = [];
    bar(exRxns) = [];
    tmpSFConsist{j}(3)=sum(foo);
    tmpSFConsist{j}(4)=sum(foo)/(sum(foo) + sum(bar));

    fluxConsistentRxnBool(exRxns) = [1];
    fluxInConsistentRxnBool(exRxns) = [0];
    
    try
        [SConsistentMetBool,SConsistentRxnBool,SInConsistentMetBool,SInConsistentRxnBool,unknownSConsistencyMetBool,unknownSConsistencyRxnBool]=...
                        findStoichConsistentSubset(model);

        foo = SConsistentRxnBool;
        bar = SInConsistentRxnBool;
        foo(exRxns) = [];
        bar(exRxns) = [];

        tmpSFConsist{j}(1)=sum(foo);
        tmpSFConsist{j}(2)=sum(foo)/(sum(foo) + sum(bar));

        SConsistentRxnBool(exRxns) = [1];
        SInConsistentRxnBool(exRxns) = [0];
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% TAKE TWO: REMOVING FLUX INCONSISTENT
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    model = modelOri;
    
    [grRatio, grRateKO, grRateWT, hasEffect, delRxn, fluxSolution] = singleRxnDeletion(model);
    
    toRemove = model.rxns(fluxInConsistentRxnBool);
    essential = model.rxns(hasEffect);
    toRemove_unessential = setdiff(toRemove,essential);
    
    model = removeRxns(model,toRemove_unessential,irrevFlag,metFlag);
    [model, unusedExchanges] = findUnusedExchangeReactions(model);
    
    [BlockedRxns] = identifyFastBlockedRxns(model);
    model = removeRxns(model,setdiff(BlockedRxns,essential),irrevFlag,metFlag);
    
    model = removeUnusedGenes(model);
    bioIndex = findRxnIDs(model, model.biomassRxnAbbr);
    model.rxns{bioIndex} = 'biomass';
    model.biomassRxnAbbr = 'biomass';
    
    resultsFileName=[resultsPath filesep model.modelID];
    save([resultsFileName '_consistent.mat'],'model')

    %fn = [resultsPath filesep model.modelID '_consistent.mat'];
    %load(fn);

    tmpStats{j}(4)=length(model.rxns);
    tmpStats{j}(5)=length(model.mets);
    tmpStats{j}(6)=length(model.genes);

    [atpFluxAerobic, atpFluxAnaerobic] = testATP(model);
    tmpATP{j}(3)=atpFluxAerobic;
    tmpATP{j}(4)=atpFluxAnaerobic;
    
    model = findSExRxnInd(model,[],printLevel-1);
    exRxns = find(~model.SIntRxnBool);
    model.lb(exRxns,1)=-1000;
    model.ub(exRxns,1)=1000;
    % remove forced flux constraints
    model.lb(model.lb>0)=0;
    % remove bounds on ATP production
    model=changeRxnBounds(model,'DM_atp_c_',1000,'u');

    % test stochiometric and flux consistency
    [fluxConsistentMetBool, fluxConsistentRxnBool, fluxInConsistentMetBool, fluxInConsistentRxnBool] = findFluxConsistentSubset(model,param);
    % exclude exchange and demand reactions
    foo = fluxConsistentRxnBool;
    bar = fluxInConsistentRxnBool;
    foo(exRxns) = [];
    bar(exRxns) = [];
    tmpSFConsist{j}(7)=sum(foo);
    tmpSFConsist{j}(8)=sum(foo)/(sum(foo) + sum(bar));

    fluxConsistentRxnBool(exRxns) = [1];
    fluxInConsistentRxnBool(exRxns) = [0];

    try
        [SConsistentMetBool,SConsistentRxnBool,SInConsistentMetBool,SInConsistentRxnBool,unknownSConsistencyMetBool,unknownSConsistencyRxnBool]=...
                        findStoichConsistentSubset(model);

        SConsistentRxnBool(exRxns) = [];
        SInConsistentRxnBool(exRxns) = [];
        tmpSFConsist{j}(5)=sum(SConsistentRxnBool);
        tmpSFConsist{j}(6)=sum(SConsistentRxnBool)/(sum(SConsistentRxnBool) + sum(SInConsistentRxnBool));
    end

    stats(j,1)=tmpStats{j}(1);
    stats(j,2)=tmpStats{j}(2);
    stats(j,3)=tmpStats{j}(3);
    atp(j,1)=tmpATP{j}(1);
    atp(j,2)=tmpATP{j}(2);
    SFconsist(j,1)=tmpSFConsist{j}(1);
    SFconsist(j,2)=tmpSFConsist{j}(2);
    if SFconsist(j,2)==0
        % assume stoichiometric consistency has to be at least as high as flux consistency
        SFconsist(j,1)=tmpSFConsist{j}(3);
        SFconsist(j,2)=tmpSFConsist{j}(4);
    end
    SFconsist(j,3)=tmpSFConsist{j}(3);
    SFconsist(j,4)=tmpSFConsist{j}(4);



    stats(j,4)=tmpStats{j}(4);
    stats(j,5)=tmpStats{j}(5);
    stats(j,6)=tmpStats{j}(6);
    atp(j,3)=tmpATP{j}(3);
    atp(j,4)=tmpATP{j}(4);
    SFconsist(j,5)=tmpSFConsist{j}(5);
    SFconsist(j,6)=tmpSFConsist{j}(6);
    if SFconsist(j,6)==0
        % assume stoichiometric consistency has to be at least as high as flux consistency
        SFconsist(j,5)=tmpSFConsist{j}(7);
        SFconsist(j,6)=tmpSFConsist{j}(8);
    end
    SFconsist(j,7)=tmpSFConsist{j}(7);
    SFconsist(j,8)=tmpSFConsist{j}(8);
end

save([resultsPath filesep 'Model_statistics'],'stats','-v7.3');
save([resultsPath filesep 'ATP_production'],'atp','-v7.3');
save([resultsPath filesep 'Stoch_Flux_Consistency'],'SFconsist','-v7.3');


dInfo = fileList(~[fileList.isdir]);
models={dInfo.name};
models=models';

table={'Model_ID','Stoichiometric consistency','Flux consistency', 'Stoichiometric consistency post-removal','Flux consistency post-removal'};
table(2:length(models)+1,1)=models;
table(2:length(models)+1,2)=cellstr(num2str(SFconsist(:,2)));
table(2:length(models)+1,3)=cellstr(num2str(SFconsist(:,4)));
table(2:length(models)+1,4)=cellstr(num2str(SFconsist(:,6)));
table(2:length(models)+1,5)=cellstr(num2str(SFconsist(:,8)));

cell2csv([resultsPath filesep 'Stoich_Flux_Consistency.csv'],table)

