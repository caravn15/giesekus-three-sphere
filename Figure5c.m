% Script to plot Figure 5c in the manuscript

clear all; close all; clc;
[p,n,e] = fileparts(mfilename('fullpath'));
rmpath([p filesep 'giesekus-verification-funcs'],'-end');
addpath([p filesep 'giesekus-three-sphere-funcs'],'-end');

%% setup

% get data
try
    load('data/All_Sims_Loaded_Data.mat');
catch
    RunSimulationsFigures5to8();
    LoadAndSortDataFigures5to8();
    load('data/All_Sims_Loaded_Data.mat');
end

% nbeats
nbeats = 125;

%% progress plots

tvec = 0:0.01:nbeats;
clr3 = [0,0,255]./255;
clr1 = [0,157,115]./255;
clr2 = [147,0,211]./255;

nplot = 15;
f2 = figure;
hold on;
plot(tvec(end-nplot*100:end),Xt_s2(end-nplot*100:end,1)-Xt_s2(1,1),...
    '-square','Color',clr3,'LineWidth',1.5,'MarkerIndices',...
    10:10:length(tvec(end-nplot*100:end)),'MarkerFaceColor',clr3)
plot(tvec(end-nplot*100:end),Xt_s2(end-nplot*100:end,612)-Xt_s2(1,1),...
    '-o','Color',clr1,'LineWidth',1.5,'MarkerIndices',...
    5:4:length(tvec(end-nplot*100:end)),'MarkerSize',4,'MarkerFaceColor',clr1)
plot(tvec(end-nplot*100:end),Xt_s2(end-nplot*100:end,end)-Xt_s2(1,1),...
    '-x','Color',clr2,'LineWidth',1.5,'MarkerIndices',...
    5:20:length(tvec(end-nplot*100:end)),'MarkerSize',10,'MarkerFaceColor',clr2)
plot(tvec(end-nplot*100:end),Xt_s2(end-nplot*100:end,end)-Xt_s2(1,1),...
    'Color',clr2,'LineWidth',1.5)
xlabel('$t$','Interpreter','latex');
ylabel('$X^{[2]}_1(t)-X^{[2]}_1(0)$','Interpreter','latex');
xlim([nbeats-nplot,nbeats])
ylim([0.23,0.41])
box on;
set(gca,'FontSize',15);
ld = legend({'Wi = $0$, $\alpha = 0$','Wi = $5$, $\alpha = 0.25$','Wi = $5$, $\alpha = 0.5$'},'Interpreter','latex');
ld.Location = 'southeast';
f2.Position = [188,131,1122,272];
ld.Box = 0;
