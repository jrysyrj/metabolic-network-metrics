clc
clear

basePath='/home/jrys/orcd/pool/metabolic-network-metrics/planets'
resultsPath=[basePath '/results/networks_MAT_clean']
mkdir(resultsPath)
printLevel=0;

fileList = dir([basePath '/results/networks_MAT/*.mat'])
fileNames = cell(length(fileList),1);
for k = 1:length(fileList) 
    fileNames{k} = strcat(fileList(k).folder,'/',fileList(k).name);
end 
fileNames = sort(fileNames);

allMets = {};
modelIDs = []; 

for k = 1:length(fileList)
    fn = fileNames{k}
    load(fn);

    modelIDs = [modelIDs, {model.modelID}];

    m = sort(model.mets);
    for i = 1:length(m) 
        if ~any(strcmp(allMets,m{i})) 
            allMets = [allMets, m{i}];
        end 
    end 
end 

allMets = sort(allMets);

for k = 1:length(fileList)
    
    fn = fileNames{k}
    load(fn);

    [rows,cols,vals] = find(model.S);
    newrows = zeros(size(rows));
    for i = 1:length(rows)
        m = model.mets(rows(i));
        newi = find(strcmp(allMets,m));
        newrows(i) = newi(1);
    end
    model.S = sparse(newrows,cols,vals,length(allMets),size(model.S,2));
    model.mets = allMets;

    [~, nRxn] = size(model.S);
    
    removedRxnInd = [];
    keptRxnInd = [];
    oneToN = 1:nRxn;
    
    % vanilla forward and reverse half stoichiometric matrices
    F        = - model.S;
    F(F < 0) = 0;
    R        = model.S;
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
                [i, ~,v] = find(model.S(:,keptOneRxnInd));
                mets = model.mets(i);
                for k = 1:length(i)
                    fprintf("(%d) %s", v(k) ,mets{k})
                    if k ~= length(i)
                        fprintf(" + ")
                    end
                end
                fprintf("\n")
                
                for k = 1:length(removedRxnInds)
                    removedOneRxnInd = removedRxnInds(k);
                    fprintf('%s (%d) \t', 'Duplicate: ',removedOneRxnInd);
                    [i, ~,v] = find(model.S(:,removedOneRxnInd));
                    mets = model.mets(i);
                    for k = 1:length(i)
                        fprintf("(%d) %s", v(k) ,mets{k})
                        if k ~= length(i)
                            fprintf(" + ")
                        end
                    end
                    fprintf("\n")
                end
    
                removedRxnInd = [removedRxnInd; removedRxnInds'];
                keptRxnInd = [keptRxnInd; keptOneRxnInd];
            end
        end
    end
    
    model.S(:,removedRxnInd) = [];
    workspace.model = model;
    save([resultsPath filesep model.modelID '.mat'],'-struct','workspace')
end
