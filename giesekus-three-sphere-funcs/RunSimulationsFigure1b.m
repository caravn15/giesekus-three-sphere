function [Xt,dt,nbeats] = RunSimulationsFigure1b()

close all; clear all; clc;

% parameters
chi    = pi/2;
Amp    = 0.8;
ep     = 0.01;
Wi     = 0;
alph   = 0;
dt     = 0.01;
tmin   = 0;
kappa  = 0.8;
nbeats = 1;
lr     = 0.2;
ng     = 14;

% swimmer initial condition
X0(:,1) = [1-lr;1];
X0(:,2) = [1;1];
X0(:,3) = [1+lr;1];


%% run simulation

[~,Xt,~,~,~,~,~,~,~,~,~,~] = RunGiesekusThreeSphereSwimmer(...
    X0,kappa,Amp,chi,lr,ep,alph,Wi,dt,nbeats,tmin,ng);

%% outputs

% save data
save('data/data_fig_1b.mat');

end