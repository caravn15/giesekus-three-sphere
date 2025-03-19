% Script to plot Figure 6 in the manuscript

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

plotStreams = 1;
refFactor   = 3;
Ulim        = [];

f3 = figure;
subplot(3,2,1);
load('data/set_1_sim_0.mat');
Utime_ind   = 1;
PlotVelocity(Utime_ind,Xt,ft,Utime,refFactor,plotStreams,Ulim);

subplot(3,2,2);
load('data/set_5_sim_0.mat');
Utime_ind   = 101;
PlotVelocity(Utime_ind,Xt,ft,Utime,refFactor,plotStreams,Ulim);

subplot(3,2,3);
load('data/set_1_sim_612.mat');
Utime_ind   = 1;
PlotVelocity(Utime_ind,Xt,ft,Utime,refFactor,plotStreams,Ulim);

subplot(3,2,4);
load('data/set_5_sim_612.mat');
Utime_ind   = 101;
PlotVelocity(Utime_ind,Xt,ft,Utime,refFactor,plotStreams,Ulim);

subplot(3,2,5);
load('data/set_1_sim_624.mat');
Utime_ind   = 1;
PlotVelocity(Utime_ind,Xt,ft,Utime,refFactor,plotStreams,Ulim);

subplot(3,2,6);
load('data/set_5_sim_624.mat');
Utime_ind   = 101;
PlotVelocity(Utime_ind,Xt,ft,Utime,refFactor,plotStreams,Ulim);

f3.Position = [0,0,612,911];

f4 = figure;
subplot(3,2,1);
load('data/set_1_sim_0.mat');
Utime_ind   = 1;
PlotVelocityZoomed(Utime_ind,Xt,ft,Utime,refFactor,plotStreams,Ulim);

subplot(3,2,2);
load('data/set_5_sim_0.mat');
Utime_ind   = 101;
PlotVelocityZoomed(Utime_ind,Xt,ft,Utime,refFactor,plotStreams,Ulim);

subplot(3,2,3);
load('data/set_1_sim_612.mat');
Utime_ind   = 1;
PlotVelocityZoomed(Utime_ind,Xt,ft,Utime,refFactor,plotStreams,Ulim);

subplot(3,2,4);
load('data/set_5_sim_612.mat');
Utime_ind   = 101;
PlotVelocityZoomed(Utime_ind,Xt,ft,Utime,refFactor,plotStreams,Ulim);

subplot(3,2,5);
load('data/set_1_sim_624.mat');
Utime_ind   = 1;
PlotVelocityZoomed(Utime_ind,Xt,ft,Utime,refFactor,plotStreams,Ulim);

subplot(3,2,6);
load('data/set_5_sim_624.mat');
Utime_ind   = 101;
PlotVelocityZoomed(Utime_ind,Xt,ft,Utime,refFactor,plotStreams,Ulim);

f4.Position= [0,0,486,719];