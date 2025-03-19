% Script to plot Figure 8 in the manuscript

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

%% plot

refFactor   = 10;

f6 = figure(12);
clf;

subplot(3,1,1)
load('data/set_5_sim_601.mat');
Ttime_ind   = 101;
PlotViscosity(Ttime_ind,Xt,ft,Ttime,refFactor,[],alph,Wi);
pbaspect([3 1 1])

subplot(3,1,2)
load('data/set_5_sim_612.mat');
Ttime_ind   = 101;
PlotViscosity(Ttime_ind,Xt,ft,Ttime,refFactor,[],alph,Wi);
pbaspect([3 1 1])

subplot(3,1,3)
load('data/set_5_sim_624.mat');
Ttime_ind   = 101;
PlotViscosity(Ttime_ind,Xt,ft,Ttime,refFactor,[],alph,Wi);
pbaspect([3 1 1])


f6.Position = [-1864,-16,1790,939];