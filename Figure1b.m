% Script to plot Figure 1(b) in the manuscript

clear all; close all; clc;
[p,~,~] = fileparts(mfilename('fullpath'));
addpath([p '/giesekus-three-sphere-funcs']);
rmpath([p filesep 'giesekus-verification-funcs'],'-end');

%% setup

% get data
try
    load('data/data_fig_1b.mat');
catch
    [Xt,dt,nbeats] = RunSimulationsFigure1b();
end

% spacing for sphere data points
step = 5;

% define time vector
tvec = 0:dt:nbeats;

%% plot

figure;
hold on;
scatter(squeeze(Xt(1,1,1:step:101))-squeeze(Xt(1,2,1:step:101)),...
    tvec(1:step:101),'k','filled')
scatter(squeeze(Xt(1,2,1:step:101))-squeeze(Xt(1,2,1:step:101)),...
    tvec(1:step:101),'k','filled')
scatter(squeeze(Xt(1,3,1:step:101))-squeeze(Xt(1,2,1:step:101)),...
    tvec(1:step:101),'k','filled')
set(gca,'FontSize',15);
box on;
xlabel('relative position of central sphere','Interpreter','latex');
ylabel('beat','Interpreter','latex');
axis square;
