clc
clear

addpath('/home/jrys/cobratoolbox')
initCobraToolbox
[solverOK,solverInstalled] = changeCobraSolver('gurobi')

basePath='/home/jrys/orcd/pool/metabolic-network-metrics/planets'
resultsPath=[basePath '/results/consistent_planets']
mkdir(resultsPath)
printLevel=0;

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


fileList = dir([basePath '/results/networks_MAT_clean/*.mat']);
fileNames = cell(length(fileList),1);
for k = 1:length(fileList) 
    fileNames{k} = strcat(fileList(k).folder,'/',fileList(k).name);
end 
fileNames = sort(fileNames);

num_rxns = [];
modelIDs = []; 
allMets = {};
allRxns = {};

S = [];

for k = 1:length(fileList)
    fn = fileNames{k}
    load(fn);

    S = [S model.S]
    num_rxns = [num_rxns; size(model.S,2)];
    modelIDs = [modelIDs, {model.modelID}];
    allMets = model.mets;
end

[~, nRxn] = size(S);

removedRxnInd = [];
keptRxnInd = [];
oneToN = 1:nRxn;

% vanilla forward and reverse half stoichiometric matrices
F        = - S;
F(F < 0) = 0;
R        = S;
R(R < 0) = 0;

A = F + R;  % invariant to direction of reaction

% detect the cols of A that are identical upto scalar multiplication
% divide each col by the sum of each row.
sumA1 = sum(A, 1);
sumA1(sumA1 == 0) = 1;
normalA1 = A * diag(1 ./ sumA1);

[~, ia, ic] = unique(normalA1', 'rows', 'stable');
nDuplicates = length(ic) - length(ia);
if nDuplicates > 0
    fprintf('%u%s\n', nDuplicates, ' duplicate reaction(s) (up to orientation)')
    for n = 1:nRxn
        bool = (ic == n);
        if nnz(bool) > 1
            ind = oneToN(bool);

            keptOneRxnInd = ind(1);
            removedRxnInds = ind(2:end);

            if length(ind) > 2
                warning(['Reaction: ' num2str(ind(1)) ' has more than one replicate'])
            end

            fprintf('%s (%d) \t', '     Keep: ',keptOneRxnInd);
            [i, ~,v] = find(S(:,keptOneRxnInd));
            mets = model.mets(i);
            rxnstr = [];
            for k = 1:length(i)
                fprintf('(%d) %s', v(k) ,mets{k});
                rxnstr = [rxnstr '(' num2str(v(k)) ')_' mets{k}];
                if k ~= length(i)
                    fprintf(" + ");
                    rxnstr = [rxnstr '_+_'];
                end
            end
            fprintf("\n")

            if length(mets) == 1
                rxnstr = ['EX_' rxnstr];
            end

            if ~any(strcmp(allRxns,rxnstr)) 
                allRxns = [allRxns, rxnstr];
            end 
            
            for k = 1:length(removedRxnInds)
                removedOneRxnInd = removedRxnInds(k);
                fprintf('%s (%d) \t', 'Duplicate: ',removedOneRxnInd);
                [i, ~,v] = find(S(:,removedOneRxnInd));
                mets = model.mets(i);
                for k = 1:length(i)
                    fprintf("(%d) %s", v(k) ,mets{k})
                    if k ~= length(i)
                        fprintf(" + ")
                    end
                end
                fprintf("\n")

                S(:,removedOneRxnInd) = S(:,keptOneRxnInd);
            end
        end
    end
end

for j = 1:size(S,2)
    if sum(sign(S(:,j)) == 1,"all")/nnz(S(:,j)) == 1
        j
        S(:,j) = -1*S(:,j);
    end
end

model_rxns = {};

for n = 1:nRxn
    [i, ~,v] = find(S(:,n));
    mets = model.mets(i);
    
    
    
    rxnstr = [];
    for k = 1:length(i)
        rxnstr = [rxnstr '(' num2str(v(k)) ')_' mets{k}];
        if k ~= length(i)
            rxnstr = [rxnstr '_+_'];
        end
    end

    if length(mets) == 1
        rxnstr = ['EX_' mets{1}];
    end
    
    model_rxns = [model_rxns, rxnstr];
    if ~any(strcmp(allRxns,rxnstr)) 
        allRxns = [allRxns, rxnstr];
    end 

end

allRxns = sort(allRxns);

cs = cumsum([1; num_rxns]);

tmpStats={};
tmpSFConsist={};

FINAL_STOICH = zeros(length(allMets),length(allRxns),length(modelIDs));

for i = 1:length(cs)-1
    tempS = S(:,cs(i):cs(i+1)-1);
    temprxns = model_rxns(cs(i):cs(i+1)-1);

    [rows,cols,vals] = find(tempS);
    newcols = zeros(size(cols));
    for j = 1:length(cols)
        r = temprxns(cols(j));
        newj = find(strcmp(allRxns,r));
        newcols(j) = newj(1);
    end
    tempS = sparse(rows,newcols,vals,length(allMets),length(allRxns));

    model = struct();

    [rows,cols,vals] = find(tempS);
    rows = unique(rows);
    cols = unique(cols);
    model.S         = tempS(rows,cols);
    model.mets      = allMets(rows)';
    %model.metFormulas = model.mets;
    model.rxns      = allRxns(cols)';
    model.ub        =  1000*ones(length(model.rxns),1);
    model.lb        = -1000*ones(length(model.rxns),1);
    model.modelID   = modelIDs{i};
    modelIDs{i}
    model.csense    = repmat('E',length(model.mets),1);
    model.b         = zeros(length(model.mets),1);
    model.c         = zeros(length(model.rxns),1);
    model.osenseStr = 'min';
    model.genes = {};
    
    tmpStats{i}(1)=length(model.rxns);
    tmpStats{i}(2)=length(model.mets);
    
    modelOri = model;
    model = findSExRxnInd(model,[],printLevel-1);
    exRxns = find(~model.SIntRxnBool);
    
    % remove constraints on exchange reactions
    model.lb(exRxns,1)=-1000;
    model.ub(exRxns,1)=1000;
    
    % remove forced flux constraints
    model.lb(model.lb>0)=0;

    % test stochiometric and flux consistency
    [fluxConsistentMetBool, fluxConsistentRxnBool, fluxInConsistentMetBool, fluxInConsistentRxnBool] = findFluxConsistentSubset(model,param);
    % exclude exchange and demand reactions
    model = findSExRxnInd(model,[],printLevel-1);
    exRxns = find(~model.SIntRxnBool);
    
    foo = fluxConsistentRxnBool;
    bar = fluxInConsistentRxnBool;
    foo(exRxns) = [];
    bar(exRxns) = [];
    tmpSFConsist{i}(3)=sum(foo);
    tmpSFConsist{i}(4)=sum(foo)/(sum(foo) + sum(bar));

    fluxConsistentRxnBool(exRxns) = [1];
    fluxInConsistentRxnBool(exRxns) = [0];
    
    try
        [SConsistentMetBool,SConsistentRxnBool,SInConsistentMetBool,SInConsistentRxnBool,unknownSConsistencyMetBool,unknownSConsistencyRxnBool]=...
                        findStoichConsistentSubset(model);

        foo = SConsistentRxnBool;
        bar = SInConsistentRxnBool;
        foo(exRxns) = [];
        bar(exRxns) = [];

        tmpSFConsist{i}(1)=sum(foo);
        tmpSFConsist{i}(2)=sum(foo)/(sum(foo) + sum(bar));

        SConsistentRxnBool(exRxns) = [1];
        SInConsistentRxnBool(exRxns) = [0];
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% TAKE TWO: REMOVING FLUX INCONSISTENT
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    model = modelOri;
    
    model = removeRxns(model,model.rxns(fluxInConsistentRxnBool),irrevFlag,metFlag);
    %[model, unusedExchanges] = findUnusedExchangeReactions(model);
    
    [BlockedRxns] = identifyFastBlockedRxns(model);
    model = removeRxns(model,BlockedRxns,irrevFlag,metFlag);
    
    resultsFileName=[resultsPath filesep model.modelID];
    save([resultsFileName '_consistent.mat'],'model')

    tmpStats{i}(3)=length(model.rxns);
    tmpStats{i}(4)=length(model.mets);

    model = findSExRxnInd(model,[],printLevel-1);
    exRxns = find(~model.SIntRxnBool);
    model.lb(exRxns,1)=-1000;
    model.ub(exRxns,1)=1000;
    % remove forced flux constraints
    model.lb(model.lb>0)=0;

    % test stochiometric and flux consistency
    [fluxConsistentMetBool, fluxConsistentRxnBool, fluxInConsistentMetBool, fluxInConsistentRxnBool] = findFluxConsistentSubset(model,param);
    % exclude exchange and demand reactions
    foo = fluxConsistentRxnBool;
    bar = fluxInConsistentRxnBool;
    foo(exRxns) = [];
    bar(exRxns) = [];
    tmpSFConsist{i}(7)=sum(foo);
    tmpSFConsist{i}(8)=sum(foo)/(sum(foo) + sum(bar));

    fluxConsistentRxnBool(exRxns) = [1];
    fluxInConsistentRxnBool(exRxns) = [0];

    try
        [SConsistentMetBool,SConsistentRxnBool,SInConsistentMetBool,SInConsistentRxnBool,unknownSConsistencyMetBool,unknownSConsistencyRxnBool]=...
                        findStoichConsistentSubset(model);

        SConsistentRxnBool(exRxns) = [];
        SInConsistentRxnBool(exRxns) = [];
        tmpSFConsist{i}(5)=sum(SConsistentRxnBool);
        tmpSFConsist{i}(6)=sum(SConsistentRxnBool)/(sum(SConsistentRxnBool) + sum(SInConsistentRxnBool));
    end

    stats(i,1)=tmpStats{i}(1);
    stats(i,2)=tmpStats{i}(2);
    SFconsist(i,1)=tmpSFConsist{i}(1);
    SFconsist(i,2)=tmpSFConsist{i}(2);
    if SFconsist(i,2)==0
        % assume stoichiometric consistency has to be at least as high as flux consistency
        SFconsist(i,1)=tmpSFConsist{i}(3);
        SFconsist(i,2)=tmpSFConsist{i}(4);
    end
    SFconsist(i,3)=tmpSFConsist{i}(3);
    SFconsist(i,4)=tmpSFConsist{i}(4);


    stats(i,3)=tmpStats{i}(3);
    stats(i,4)=tmpStats{i}(4);
    SFconsist(i,5)=tmpSFConsist{i}(5);
    SFconsist(i,6)=tmpSFConsist{i}(6);
    if SFconsist(i,6)==0
        % assume stoichiometric consistency has to be at least as high as flux consistency
        SFconsist(i,5)=tmpSFConsist{i}(7);
        SFconsist(i,6)=tmpSFConsist{i}(8);
    end
    SFconsist(i,7)=tmpSFConsist{i}(7);
    SFconsist(i,8)=tmpSFConsist{i}(8);
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

