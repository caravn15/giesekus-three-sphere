% Script to plot Figure 2(b) in the manuscript

clear all; close all; clc;
[p,~,~] = fileparts(mfilename('fullpath'));
addpath([p '/giesekus-three-sphere-funcs']);
rmpath([p filesep 'giesekus-verification'],'-end');

%% setup

% parameters
ep    = 0.01;
lx    = [0,2];
ly    = [0,2];
ng    = 14;

% get data
try
    load('data/data_fig_2b.mat');
catch
    [vels] = RunSimulationsFigure2b();
end

%% calculations

% recreate mesh
xv      = linspace(lx(1),lx(2),2*ng+1);
yv      = linspace(ly(1),ly(2),2*ng+1);
[xv,yv] = meshgrid(xv,yv);

%% plot
figure;
subplot(1,3,1);
[u,v] = ExtractVectorComponents(vels{1}.uN);
quiver(xv,yv,reshape(u,size(xv)),reshape(v,size(xv)),10)
xlim([0,2])
ylim([0,2])
axis square;
box on;
set(gca,'XTick',[])
set(gca,'YTick',[])

subplot(1,3,2);
[u,v] = ExtractVectorComponents(vels{1}.U);
quiver(xv,yv,reshape(u,size(xv)),reshape(v,size(xv)),1)
xlim([0,2])
ylim([0,2])
axis square;
box on;
set(gca,'XTick',[])
set(gca,'YTick',[])

subplot(1,3,3);
[u,v] = ExtractVectorComponents(vels{1}.u);
quiver(xv,yv,reshape(u,size(xv)),reshape(v,size(xv)),10)
xlim([0,2])
ylim([0,2])
axis square;
box on;
set(gca,'XTick',[])
set(gca,'YTick',[])