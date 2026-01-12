%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ECOLI ALL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear 
clc

load('/home/jrys/orcd/pool/metabolic-network-metrics/bacteria/results/MetRxnInfo.mat')

basePath='/home/jrys/orcd/pool/metabolic-network-metrics/bacteria/results/';
resultsPath=[basePath 'parsed_ecoli'];
fileList = dir([resultsPath filesep '*.txt']);
fileNames = cell(length(fileList),1);
for k = 1:length(fileList) 
    fileNames{k} = strcat(fileList(k).folder,'/',fileList(k).name);
end 
fileNames = sort(fileNames);

data = zeros(length(fileNames),length(rxns));
data2 = zeros(length(fileNames),length(mets));

for k = 1:length(fileList)
    fn = fileNames{k}
    M = readmatrix(fn);
    r = sort(unique(M(:,2)));
    m = sort(unique(M(:,1)));
    for j = 1:length(r) 
        data(k,r(j)) = 1;
    end
    for j = 1:length(m) 
        data2(k,m(j)) = 1;
    end
end 

fn_out = '/home/jrys/orcd/pool/metabolic-network-metrics/bacteria/results/ecoli_rxns_present.csv';
dlmwrite(fn_out,data,'delimiter',',');

fn_out = '/home/jrys/orcd/pool/metabolic-network-metrics/bacteria/results/ecoli_mets_present.csv';
dlmwrite(fn_out,data2,'delimiter',',');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ECOLI INTERNAL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc

load('/home/jrys/orcd/pool/metabolic-network-metrics/bacteria/results/MetRxnInfo.mat')

basePath='/home/jrys/orcd/pool/metabolic-network-metrics/bacteria/results/';
resultsPath=[basePath 'parsed_ecoli_int'];
fileList = dir([resultsPath filesep '*.txt']);
fileNames = cell(length(fileList),1);
for k = 1:length(fileList) 
    fileNames{k} = strcat(fileList(k).folder,'/',fileList(k).name);
end 
fileNames = sort(fileNames);

data = zeros(length(fileNames),length(rxns));
data2 = zeros(length(fileNames),length(mets));

for k = 1:length(fileList)
    fn = fileNames{k}
    M = readmatrix(fn);
    r = sort(unique(M(:,2)));
    m = sort(unique(M(:,1)));
    for j = 1:length(r) 
        data(k,r(j)) = 1;
    end
    for j = 1:length(m) 
        data2(k,m(j)) = 1;
    end
end 

fn_out = '/home/jrys/orcd/pool/metabolic-network-metrics/bacteria/results/ecoli_rxns_present_int.csv';
dlmwrite(fn_out,data,'delimiter',',');

fn_out = '/home/jrys/orcd/pool/metabolic-network-metrics/bacteria/results/ecoli_mets_present_int.csv';
dlmwrite(fn_out,data2,'delimiter',',');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SUBAGORA2 ALL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear 
clc

load('/home/jrys/orcd/pool/metabolic-network-metrics/bacteria/results/MetRxnInfo.mat')

basePath='/home/jrys/orcd/pool/metabolic-network-metrics/bacteria/results/';
resultsPath=[basePath 'parsed_subAGORA2'];
fileList = dir([resultsPath filesep '*.txt']);
fileNames = cell(length(fileList),1);
for k = 1:length(fileList) 
    fileNames{k} = strcat(fileList(k).folder,'/',fileList(k).name);
end 
fileNames = sort(fileNames);

data = zeros(length(fileNames),length(rxns));
data2 = zeros(length(fileNames),length(mets));

for k = 1:length(fileList)
    fn = fileNames{k}
    M = readmatrix(fn);
    r = sort(unique(M(:,2)));
    m = sort(unique(M(:,1)));
    for j = 1:length(r) 
        data(k,r(j)) = 1;
    end
    for j = 1:length(m) 
        data2(k,m(j)) = 1;
    end
end 

fn_out = '/home/jrys/orcd/pool/metabolic-network-metrics/bacteria/results/subAGORA2_rxns_present.csv';
dlmwrite(fn_out,data,'delimiter',',');

fn_out = '/home/jrys/orcd/pool/metabolic-network-metrics/bacteria/results/subAGORA2_mets_present.csv';
dlmwrite(fn_out,data2,'delimiter',',');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SUBAGORA2 INTERNAL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc

load('/home/jrys/orcd/pool/metabolic-network-metrics/bacteria/results/MetRxnInfo.mat')

basePath='/home/jrys/orcd/pool/metabolic-network-metrics/bacteria/results/';
resultsPath=[basePath 'parsed_subAGORA2_int'];
fileList = dir([resultsPath filesep '*.txt']);
fileNames = cell(length(fileList),1);
for k = 1:length(fileList) 
    fileNames{k} = strcat(fileList(k).folder,'/',fileList(k).name);
end 
fileNames = sort(fileNames);

data = zeros(length(fileNames),length(rxns));
data2 = zeros(length(fileNames),length(mets));

for k = 1:length(fileList)
    fn = fileNames{k}
    M = readmatrix(fn);
    r = sort(unique(M(:,2)));
    m = sort(unique(M(:,1)));
    for j = 1:length(r) 
        data(k,r(j)) = 1;
    end
    for j = 1:length(m) 
        data2(k,m(j)) = 1;
    end
end 

fn_out = '/home/jrys/orcd/pool/metabolic-network-metrics/bacteria/results/subAGORA2_rxns_present_int.csv';
dlmwrite(fn_out,data,'delimiter',',');

fn_out = '/home/jrys/orcd/pool/metabolic-network-metrics/bacteria/results/subAGORA2_mets_present_int.csv';
dlmwrite(fn_out,data2,'delimiter',',');