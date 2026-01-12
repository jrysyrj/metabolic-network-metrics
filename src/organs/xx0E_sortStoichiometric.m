clc
clear

addpath('/home/jrys/cobratoolbox')
initCobraToolbox
[solverOK,solverInstalled] = changeCobraSolver('gurobi')

basePath='/home/jrys/orcd/pool/metabolic-network-metrics/organs/results/'
resultPath    =[basePath 'parsed_organs/']
resultPath_int=[basePath 'parsed_organs_int/']
mkdir(resultPath)
mkdir(resultPath_int)
fileList = dir([basePath 'consistent_organs/*_consistent.mat']);
fileNames = cell(length(fileList),1);
for k = 1:length(fileList) 
    fileNames{k} = strcat(fileList(k).folder,'/',fileList(k).name);
end 
fileNames = sort(fileNames);

load([basePath 'MetRxnInfo.mat'])

for k = 1:length(fileNames)
    fn = fileNames{k}
    load(fn);
    
    [rows,cols,vals] = find(model.S);
    newrows = zeros(size(rows));
    newcols = zeros(size(cols));
    for i = 1:length(rows)
        m = model.mets(rows(i));
        newi = find(strcmp(mets,m));
        newrows(i) = newi(1);
    end
    for j = 1:length(cols)
        r = model.rxns(cols(j));
        newj = find(strcmp(rxns,r));
        newcols(j) = newj(1);
    end
    newS = [newrows,newcols,vals];
    newfn = strcat(resultPath,pad(string(k),3,'left','0'),'.txt');
    dlmwrite(newfn,newS,'delimiter','\t');
end 

%%
printLevel = 0;
for k = 1:length(fileNames)
    fn = fileNames{k}
    load(fn);

    model = findSExRxnInd(model,[],printLevel-1);
    exRxns = find(~model.SIntRxnBool);
    
    S = model.S;
    S(:,exRxns) = [];
    modelrxns = model.rxns;
    modelrxns(exRxns) = [];
    
    [rows,cols,vals] = find(S);
    newrows = zeros(size(rows));
    newcols = zeros(size(cols));
    for i = 1:length(rows)
        m = model.mets(rows(i));
        newi = find(strcmp(mets,m));
        newrows(i) = newi(1);
    end
    for j = 1:length(cols)
        r = modelrxns(cols(j));
        newj = find(strcmp(rxns,r));
        newcols(j) = newj(1);
    end
    newS = [newrows,newcols,vals];
    newfn = strcat(resultPath_int,pad(string(k),3,'left','0'),'.txt');
    dlmwrite(newfn,newS,'delimiter','\t');
end 
