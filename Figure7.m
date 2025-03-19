% Script to plot Figure 7 in the manuscript

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

%% plots

refFactor   = 10;

f5 = figure;
subplot(3,2,1)
load('data/set_5_sim_0.mat');
Ttime_ind   = 100;
PlotStress(Ttime_ind,Xt,ft,Ttime,refFactor,[0,0.52]);
pbaspect([3 1 1])

subplot(3,2,3)
load('data/set_5_sim_612.mat');
Ttime_ind   = 100;
PlotStress(Ttime_ind,Xt,ft,Ttime,refFactor,[0,0.52]);
pbaspect([3 1 1])

subplot(3,2,5)
load('data/set_5_sim_624.mat');
Ttime_ind   = 100;
PlotStress(Ttime_ind,Xt,ft,Ttime,refFactor,[0,0.52]);
pbaspect([3 1 1])

subplot(3,2,2)
load('data/set_5_sim_0.mat');
Ttime_ind   = 101;
PlotStress(Ttime_ind,Xt,ft,Ttime,refFactor,[0,0.62]);
pbaspect([3 1 1])

subplot(3,2,4)
load('data/set_5_sim_612.mat');
Ttime_ind   = 101;
PlotStress(Ttime_ind,Xt,ft,Ttime,refFactor,[0,0.62]);
pbaspect([3 1 1])

subplot(3,2,6)
load('data/set_5_sim_624.mat');
Ttime_ind   = 101;
PlotStress(Ttime_ind,Xt,ft,Ttime,refFactor,[0,0.62]);
pbaspect([3 1 1])

f5.Position = [-1768,18,1277,951];

