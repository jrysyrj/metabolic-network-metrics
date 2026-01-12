clc
clear

addpath('/home/jrys/cobratoolbox')
initCobraToolbox
[solverOK,solverInstalled] = changeCobraSolver('gurobi')

basePath='/home/jrys/orcd/pool/metabolic-network-metrics/bacteria'
resultsPath=[basePath '/results/']
resultsPath2=[resultsPath 'ecoliKOs/']
mkdir(resultsPath2)
printLevel=0;


fn = [resultsPath 'ECOLI_TRANSLATED.mat']

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
tmpSFConsist={};

load(fn);
modelOri = model;

geneKO = calculateGeneKOMatrix(modelOri)

[geneKOmatrix_u,ia,ic] = unique(geneKO.matrix','rows','stable');
geneKOmatrix_u = geneKOmatrix_u';
size(geneKOmatrix_u,1) == size(geneKO.matrix,1)

geneGroups = [];
for k = 1:size(geneKOmatrix_u,2)
    targetCol= geneKOmatrix_u(:,k);
    A = geneKO.matrix';
    b = targetCol';
    [~, idx] = ismember(A, b, 'rows');
    idx = find(idx);
    geneGroups = [geneGroups; join(geneKO.genes(idx),"," )];
end

writecell(geneGroups, [resultsPath2 'geneGroups.csv']); 

for j = 1:length(geneGroups)
    model = modelOri;
    
    tmpStats{j}(1)=length(model.rxns);
    tmpStats{j}(2)=length(model.mets);
    tmpStats{j}(3)=length(model.genes);

    
    toRemove = geneKO.rxns(geneKOmatrix_u(:,j));
    model = removeRxns(model,toRemove,irrevFlag,metFlag);
    model = removeUnusedGenes(model);
    
    
    modelp = model;
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
    

    toRemove = model.rxns(fluxInConsistentRxnBool);
    model = modelp;
    model = removeRxns(model,toRemove,irrevFlag,metFlag);
    model = removeUnusedGenes(model);
    
    resultsFileName = strcat(resultsPath2,'geneGroup',pad(string(j),3,'left','0'),'_consistent.mat');
    save(resultsFileName,'model')

    tmpStats{j}(4)=length(model.rxns);
    tmpStats{j}(5)=length(model.mets);
    tmpStats{j}(6)=length(model.genes);

    
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

save([resultsPath2 filesep 'Model_statistics'],'stats','-v7.3');
save([resultsPath2 filesep 'Stoch_Flux_Consistency'],'SFconsist','-v7.3');

models=1:length(geneGroups);
models = models';

table={'Model_ID','Stoichiometric consistency','Flux consistency', 'Stoichiometric consistency post-removal','Flux consistency post-removal'};
table(2:length(models)+1,1)=cellstr(num2str(models));
table(2:length(models)+1,2)=cellstr(num2str(SFconsist(:,2)));
table(2:length(models)+1,3)=cellstr(num2str(SFconsist(:,4)));
table(2:length(models)+1,4)=cellstr(num2str(SFconsist(:,6)));
table(2:length(models)+1,5)=cellstr(num2str(SFconsist(:,8)));

cell2csv([resultsPath2 filesep 'Stoich_Flux_Consistency.csv'],table)