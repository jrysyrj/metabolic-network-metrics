clc
clear

addpath('/home/jrys/cobratoolbox')
initCobraToolbox
[solverOK,solverInstalled] = changeCobraSolver('gurobi')

basePath='/home/jrys/orcd/pool/metabolic-network-metrics/organs'
resultsPath=[basePath '/results/consistent_organs']
mkdir(resultsPath)
printLevel=0;

fileList = dir([basePath '/results/organ_recons/*.mat'])
fileNames = cell(length(fileList),1);
for k = 1:length(fileList) 
    fileNames{k} = strcat(fileList(k).folder,'/',fileList(k).name);
end 
fileNames = sort(fileNames);

% param for consistency checks 
param=struct;
feasTol = getCobraSolverParams('LP', 'feasTol');
param.('feasTol')= feasTol;
param.('epsilon')=feasTol*10;
param.('modeFlag')=0;
param.('method')='fastcc';

% flags for rxn removal 
irrevFlag=0;
metFlag=1;

    
tmpStats={};
tmpSFConsist={};


for j = 1:length(fileNames)
    fn = fileNames{j}
    load(fn);
    
    [model, removedRxnInd, keptRxnInd] = checkDuplicateRxn(model);
    
    if any(strcmp(model.rxns,'UGT1A10r'))
        model = addReaction(model, 'UGT1A10r',...
            'metaboliteList', {'bilirub[r]', 'udpglcur[r]', 'bilglcur[r]', 'mudp[r]'},...
            'stoichCoeffList', [-1; -1; 1; 1]);
    end
    if any(strcmp(model.rxns,'PMEVKc'))
        model = addReaction(model, 'PMEVKc',...
            'metaboliteList', {'atp[c]', '5pmev[c]', 'adp[c]', '5dpmev[c]'},...
            'stoichCoeffList', [-1; -1; 1; 1]);
    end
    
    if any(strcmp(model.rxns,'DPMVDc'))
        model = addReaction(model, 'DPMVDc',...
            'metaboliteList', {'atp[c]', '5dpmev[c]', 'adp[c]', 'pi[c]','co2[c]','ipdp[c]'},...
            'stoichCoeffList', [-1; -1; 1; 1; 1; 1]);
    end 
    
    if any(strcmp(model.rxns,'DASCBR'))
        model = addReaction(model, 'DASCBR',...
            'metaboliteList', {'nadph[c]', 'dhdascb[c]', 'nadp[c]', 'ascb_L[c]'},...
            'stoichCoeffList', [-1; -1; 1; 1]);
    end 
    
    if any(strcmp(model.rxns,'CYSAMOe[luI]'))
        model = addReaction(model, 'CYSAMOe[luI]',...
            'metaboliteList', {'o2[luI]', 'cysam[luI]', 'h[luI]', 'hyptaur[luI]'},...
            'stoichCoeffList', [-1; -1; 1; 1]);
    end
    
    if any(strcmp(model.rxns,'HC02196c'))
        model = addReaction(model, 'HC02196c',...
            'metaboliteList', {'gly[c]', 'urscholcoa[c]', 'coa[c]', 'HC02196[c]'},...
            'stoichCoeffList', [-1; -1; 1; 1]);
    end
    
    if any(strcmp(model.rxns,'HC02195c'))
        model = addReaction(model, 'HC02195c',...
            'metaboliteList', {'taur[c]', 'urscholcoa[c]', 'coa[c]', 'HC02195[c]'},...
            'stoichCoeffList', [-1; -1; 1; 1]);
    end
    
    if any(strcmp(model.rxns,'HMR_1750'))
        model = addReaction(model, 'HMR_1750',...
            'metaboliteList', {'h[m]', 'nadph[m]', 'M00746[m]', 'h2o[m]','nadp[m]','M02977[m]'},...
            'stoichCoeffList', [-2; -1; -1; 1; 1; 1]);
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
    
    missingFormulaMetsIdx = [];
    for k = 1:length(model.mets)
        if length(model.metFormulas{k}) == 0
            missingFormulaMetsIdx = [missingFormulaMetsIdx; k];
        elseif contains(model.metFormulas{k},'R') 
            missingFormulaMetsIdx = [missingFormulaMetsIdx; k];
        end
    end
    
    [~,rxnWithUnknownFormulas_idxs,~] = find(model.S(missingFormulaMetsIdx,:));
    
    rxnsWithUnknownFormulas = model.rxns(unique(rxnWithUnknownFormulas_idxs));
    
    exclude = {'TECAAE','TECAGE','TECAUE','PROTEIN_BS','B3GNT313g','FUT17g'};
    exclude = [exclude, rxnsWithUnknownFormulas'];
    exclude_idxs = [];
    for k = 1:length(exclude)
        exclude_idxs = [exclude_idxs; find(strcmp(model.rxns,exclude{k}))];
    end
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
    
    modelOri = model;
    model = findSExRxnInd(model,[],printLevel-1);
    exRxns = find(~model.SIntRxnBool);
    exRxnsDEGRx = sort(unique([exRxns; find(endsWith(model.rxns,'DEGRx')) ]));
    
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
    exRxnsDEGRx = sort(unique([exRxns; find(endsWith(model.rxns,'DEGRx')) ]));
    
    toRemove = {};
    try
        [SConsistentMetBool,SConsistentRxnBool,SInConsistentMetBool,SInConsistentRxnBool,unknownSConsistencyMetBool,unknownSConsistencyRxnBool]=...
                        findStoichConsistentSubset(model);

        foo = SConsistentRxnBool;
        bar = SInConsistentRxnBool;
        
        %foo(unique([exRxnsDEGRx; exclude_idxs])) = [];
        %bar(unique([exRxnsDEGRx; exclude_idxs])) = [];
        
        foo(exRxnsDEGRx) = [];
        bar(exRxnsDEGRx) = [];

        tmpSFConsist{j}(1)=sum(foo);
        tmpSFConsist{j}(2)=sum(foo)/(sum(foo) + sum(bar));

        %SConsistentRxnBool(unique([exRxnsDEGRx; exclude_idxs])) = [1];
        %SInConsistentRxnBool(unique([exRxnsDEGRx; exclude_idxs])) = [0];
        
        SConsistentRxnBool(exRxnsDEGRx) = [1];
        SInConsistentRxnBool(exRxnsDEGRx) = [0];
        
        toRemove = model.rxns(SInConsistentRxnBool);
    end
    
    model = modelOri;
    model = removeRxns(model,toRemove,irrevFlag,metFlag);
    
    modelOri = model;
    model = findSExRxnInd(model,[],printLevel-1);
    exRxns = find(~model.SIntRxnBool);
    exRxnsDEGRx = sort(unique([exRxns; find(endsWith(model.rxns,'DEGRx')) ]));
    
    % remove constraints on exchange reactions
    model.lb(exRxns,1)=-1000;
    model.ub(exRxns,1)=1000;
    
    % remove forced flux constraints
    model.lb(model.lb>0)=0;

    % remove bounds on ATP production
    model=changeRxnBounds(model,'DM_atp_c_',1000,'u');
    
    % test stochiometric and flux consistency
    [fluxConsistentMetBool, fluxConsistentRxnBool, fluxInConsistentMetBool, fluxInConsistentRxnBool] = findFluxConsistentSubset(model,param);
    
    foo = fluxConsistentRxnBool;
    bar = fluxInConsistentRxnBool;
    foo(exRxnsDEGRx) = [];
    bar(exRxnsDEGRx) = [];
    tmpSFConsist{j}(3)=sum(foo);
    tmpSFConsist{j}(4)=sum(foo)/(sum(foo) + sum(bar));

    fluxConsistentRxnBool(exRxnsDEGRx) = [1];
    fluxInConsistentRxnBool(exRxnsDEGRx) = [0];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% TAKE TWO: REMOVING FLUX INCONSISTENT
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    model = modelOri;
    
    toRemove = model.rxns(fluxInConsistentRxnBool);
    model = removeRxns(model,toRemove,irrevFlag,metFlag);
    [model, unusedExchanges] = findUnusedExchangeReactions(model);
    [BlockedRxns] = identifyFastBlockedRxns (model);
    model = removeRxns(model,BlockedRxns,irrevFlag,metFlag);
    
    
    
    
    
    
    modelOri = model;
    model = findSExRxnInd(model,[],printLevel-1);
    exRxns = find(~model.SIntRxnBool);
    exRxnsDEGRx = sort(unique([exRxns; find(endsWith(model.rxns,'DEGRx')) ]));
    
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
    exRxnsDEGRx = sort(unique([exRxns; find(endsWith(model.rxns,'DEGRx')) ]));
    
    toRemove = {};
    try
        [SConsistentMetBool,SConsistentRxnBool,SInConsistentMetBool,SInConsistentRxnBool,unknownSConsistencyMetBool,unknownSConsistencyRxnBool]=...
                        findStoichConsistentSubset(model);

        foo = SConsistentRxnBool;
        bar = SInConsistentRxnBool;
        
        %foo(unique([exRxnsDEGRx; exclude_idxs])) = [];
        %bar(unique([exRxnsDEGRx; exclude_idxs])) = [];
        
        foo(exRxnsDEGRx) = [];
        bar(exRxnsDEGRx) = [];

        tmpSFConsist{j}(1)=sum(foo);
        tmpSFConsist{j}(2)=sum(foo)/(sum(foo) + sum(bar));

        %SConsistentRxnBool(unique([exRxnsDEGRx; exclude_idxs])) = [1];
        %SInConsistentRxnBool(unique([exRxnsDEGRx; exclude_idxs])) = [0];
        
        SConsistentRxnBool(exRxnsDEGRx) = [1];
        SInConsistentRxnBool(exRxnsDEGRx) = [0];
        
        toRemove = model.rxns(SInConsistentRxnBool);
    end
    
    model = modelOri;
    model = removeRxns(model,toRemove,irrevFlag,metFlag);
    
    modelOri = model;
    model = findSExRxnInd(model,[],printLevel-1);
    exRxns = find(~model.SIntRxnBool);
    exRxnsDEGRx = sort(unique([exRxns; find(endsWith(model.rxns,'DEGRx')) ]));
    
    % remove constraints on exchange reactions
    model.lb(exRxns,1)=-1000;
    model.ub(exRxns,1)=1000;
    
    % remove forced flux constraints
    model.lb(model.lb>0)=0;

    % remove bounds on ATP production
    model=changeRxnBounds(model,'DM_atp_c_',1000,'u');
    
    % test stochiometric and flux consistency
    [fluxConsistentMetBool, fluxConsistentRxnBool, fluxInConsistentMetBool, fluxInConsistentRxnBool] = findFluxConsistentSubset(model,param);
    
    foo = fluxConsistentRxnBool;
    bar = fluxInConsistentRxnBool;
    foo(exRxnsDEGRx) = [];
    bar(exRxnsDEGRx) = [];
    tmpSFConsist{j}(3)=sum(foo);
    tmpSFConsist{j}(4)=sum(foo)/(sum(foo) + sum(bar));

    fluxConsistentRxnBool(exRxnsDEGRx) = [1];
    fluxInConsistentRxnBool(exRxnsDEGRx) = [0];
    
    
    
    
    
    
    model = modelOri;
    toRemove = model.rxns(fluxInConsistentRxnBool);
    model = removeRxns(model,toRemove,irrevFlag,metFlag);
    
    model = removeUnusedGenes(model);
    
    resultsFileName=[resultsPath filesep model.modelID];
    save([resultsFileName '_consistent.mat'],'model')

    tmpStats{j}(4)=length(model.rxns);
    tmpStats{j}(5)=length(model.mets);
    tmpStats{j}(6)=length(model.genes);

    
    model = findSExRxnInd(model,[],printLevel-1);
    exRxns = find(~model.SIntRxnBool);
    exRxnsDEGRx = sort(unique([exRxns; find(endsWith(model.rxns,'DEGRx')) ]));
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
    foo(exRxnsDEGRx) = [];
    bar(exRxnsDEGRx) = [];
    tmpSFConsist{j}(7)=sum(foo);
    tmpSFConsist{j}(8)=sum(foo)/(sum(foo) + sum(bar));

    fluxConsistentRxnBool(exRxnsDEGRx) = [1];
    fluxInConsistentRxnBool(exRxnsDEGRx) = [0];

    exclude_idxs = [];
    for k = 1:length(exclude)
        exclude_idxs = [exclude_idxs; find(strcmp(model.rxns,exclude{k}))];
    end
    try
        [SConsistentMetBool,SConsistentRxnBool,SInConsistentMetBool,SInConsistentRxnBool,unknownSConsistencyMetBool,unknownSConsistencyRxnBool]=...
                        findStoichConsistentSubset(model);

        %SConsistentRxnBool(unique([exRxnsDEGRx; exclude_idxs])) = [];
        %SInConsistentRxnBool(unique([exRxnsDEGRx; exclude_idxs])) = [];
        
        SConsistentRxnBool(exRxnsDEGRx) = [1];
        SInConsistentRxnBool(exRxnsDEGRx) = [0];
        
        tmpSFConsist{j}(5)=sum(SConsistentRxnBool);
        tmpSFConsist{j}(6)=sum(SConsistentRxnBool)/(sum(SConsistentRxnBool) + sum(SInConsistentRxnBool));
    end

    stats(j,1)=tmpStats{j}(1);
    stats(j,2)=tmpStats{j}(2);
    stats(j,3)=tmpStats{j}(3);
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
