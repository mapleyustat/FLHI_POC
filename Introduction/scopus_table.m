%% Generates a table relating the number of publications to search terms
% from the SCOPUS database (www.scopus.com)

clear all;
close all;
clc;

% Variables that control column indices
Hidx = 2;
MPCidx = 3;
MPCHidx = 4;
DMCidx = 5;
DMCHidx = 6;

%% First column is the year and second column the number of publications
dataDMC = load('ScopusDMC.txt'); 
dataDMCH = load('ScopusDMCHammerstein.txt'); 
dataMPC = load('ScopusMPC.txt'); 
dataMPCH = load('ScopusMPCHammerstein.txt'); 
dataH = load('ScopusHammerstein.txt'); 

%% Remove the first year (2014) since it's incomplete (data from 27/01/2014)
dataDMC(1,:) = [];
dataMPC(1,:) = [];
dataMPCH(1,:) = [];
dataH(1,:) = [];

%% Loop through the years 2013 -> 1978 and gather data
dataAll = zeros(length([2013:-1:1978]), 6);

% Assign the years 2013 to 1978
dataAll(:,1) = [2013:-1:1978]';

% Assign Hammerstein years
[A,B] = ismember(dataAll(:,1), dataH(:,1));
dataAll(A,Hidx) = dataH(B(B~=0),2);
% special case, Hammerstein goes beyond 1978 so sum up all the years
dataAll(end,Hidx) = sum(dataH(dataH(:,1) <= 1978,2))

[A,~] = ismember(dataAll(:,1), dataDMC(:,1));
dataAll(A,DMCidx) = dataDMC(:,2);

[A,~] = ismember(dataAll(:,1), dataDMCH(:,1));
dataAll(A,DMCHidx) = dataDMCH(:,2);

[A,~] = ismember(dataAll(:,1), dataMPC(:,1));
dataAll(A,MPCidx) = dataMPC(:,2);

[A,~] = ismember(dataAll(:,1), dataMPCH(:,1));
dataAll(A,MPCHidx) = dataMPCH(:,2);

dataAll

%input.tableData = dataAll;
%input.tableRowLabels = cellstr(num2str(dataAll(:,1)));
%latexTable(input)

for k=1:size(dataAll,1)
    fprintf('%g & %g & %g & %g & %g & %g \\\\ \n', dataAll(k,:))
end
