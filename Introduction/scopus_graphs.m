%% Generates graphs relating the number of publications to search terms
% from the SCOPUS database (www.scopus.com)

clear all;
close all;
clc;

% Variables which control figure output
fontsize = 12;
ylabelText = 'No. of Publications';
xlabelText = 'Year';
%fontsize = 14;
%ylabelText = 'No. de Publicações';
%xlabelText = 'Ano';

%% First column is the year and second column the number of publications
dataDMC = load('ScopusDMC.txt'); 
dataDMCH = load('ScopusDMCHammerstein.txt'); 
dataMPC = load('ScopusMPC.txt'); 
dataMPCH = load('ScopusMPCHammerstein.txt'); 
dataH = load('ScopusHammerstein.txt'); 

%% Generates figures
% DMC
h1 = figure;
bar(dataDMC(:,1), dataDMC(:,2), 'w');
%axis([1980 2015 0 45]);
ylabel(ylabelText);
xlabel(xlabelText);

xh = get(gca, 'xlabel');
yh = get(gca, 'ylabel');
set(xh, 'fontsize', fontsize)
set(yh, 'fontsize', fontsize)
set(gca, 'fontsize', fontsize)
saveTightFigure(h1, 'ScopusDMC.fig')
print(gcf,'-dpdf','ScopusDMC.pdf')

% DMC Hammerstein
h2 = figure;
bar(dataDMCH(:,1), dataDMCH(:,2), 'w');
ylabel(ylabelText);
xlabel(xlabelText);
%axis([1995 2015 0 4]);
set(gca,'XTick',[1995 2000 2005 2010 2015])

xh = get(gca, 'xlabel');
yh = get(gca, 'ylabel');
set(xh, 'fontsize', fontsize)
set(yh, 'fontsize', fontsize)
set(gca, 'fontsize', fontsize)
saveTightFigure(h2, 'ScopusDMCH.fig')
print(gcf,'-dpdf','ScopusDMCH.pdf')

% MPC
h3 = figure;
bar(dataMPC(:,1), dataMPC(:,2), 'w');
ylabel(ylabelText);
xlabel(xlabelText);

xh = get(gca, 'xlabel');
yh = get(gca, 'ylabel');
set(xh, 'fontsize', fontsize)
set(yh, 'fontsize', fontsize)
set(gca, 'fontsize', fontsize)
saveTightFigure(h3, 'ScopusMPC.fig')
print(gcf,'-dpdf','ScopusMPC.pdf')

% MPC Hammerstein
h4 = figure;
bar(dataMPCH(:,1), dataMPCH(:,2), 'w');
ylabel(ylabelText);
xlabel(xlabelText);

xh = get(gca, 'xlabel');
yh = get(gca, 'ylabel');
set(xh, 'fontsize', fontsize)
set(yh, 'fontsize', fontsize)
set(gca, 'fontsize', fontsize)
saveTightFigure(h4, 'ScopusMPCH.fig')
print(gcf,'-dpdf','ScopusMPCH.pdf')

% Hammerstein
h5 = figure;
bar(dataH(:,1), dataH(:,2), 'w');
ylabel(ylabelText);
xlabel(xlabelText);
%axis([1950 2015 0 250]);

xh = get(gca, 'xlabel');
yh = get(gca, 'ylabel');
set(xh, 'fontsize', fontsize)
set(yh, 'fontsize', fontsize)
set(gca, 'fontsize', fontsize)
saveTightFigure(h5, 'ScopusH.fig')
print(gcf,'-dpdf','ScopusH.pdf')