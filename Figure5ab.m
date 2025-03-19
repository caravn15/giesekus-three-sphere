% Script to plot Figure 5a-b in the manuscript

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

% define alph and Wi
alp_v        = linspace(0,0.5,25);
Wi_v         = linspace(0,5,25);
[Wi_M,alp_M] = meshgrid(Wi_v,alp_v);
alp_V        = alp_M(:);
Wi_V         = Wi_M(:);

% reshape data for pcolor
Ubar_M   = reshape(Ubar(:,5),size(Wi_M));
Ubar_max = max(Ubar_M(:));

Wbar_M   = reshape(Wbar,size(Wi_M));
Wbar_max = max(Wbar_M(:));

eff_M   = reshape(eff,size(Wi_M));
eff_max = max(eff_M(:));

%% plots
f1 = figure;

% Ubar
subplot(1,2,1);
hold on;
pc                   = pcolor(Wi_M(2:end,:),alp_M(2:end,:),Ubar_M(2:end,:)./Ubar_max);
pc.EdgeAlpha         = 0.2;
cb                   = colorbar;
cb.Label.Interpreter = 'latex';
cb.Label.String      = '$\bar{U}/\bar{U}_{max}$';
cb.Label.FontSize    = 14;
cb.Location          = 'northoutside';
cb.Ticks             = [0.94,0.97,1];
scatter(Wi_v(end)-(Wi_v(end)-Wi_v(end-1))/2,alp_v(end)-...
    (alp_v(end)-alp_v(end-1))/2,65,'rx');
set(gca,'FontSize',15);
axis square;
box on;
ax = gca;
ax.Layer = 'top';
xlabel('Wi','Interpreter','latex');
ylabel('$\alpha$','Interpreter','latex');
clim([0.94,1])
xlim([Wi_v(2),Wi_v(end)])
ylim([alp_v(2),alp_v(end)])
xticks([0,1,2,3,4,5])
yticks([0,0.1,0.2,0.3,0.4,0.5])

% eff
subplot(1,2,2);
hold on;
pc                   = pcolor(Wi_M(2:end,:),alp_M(2:end,:),eff_M(2:end,:)./eff_max);
pc.EdgeAlpha         = 0.2;
cb                   = colorbar;
cb.Label.Interpreter = 'latex';
cb.Label.String      = '$\mathcal{E}/\mathcal{E}_{max}$';
cb.Label.FontSize    = 14;
cb.Location          = 'northoutside';
cb.Ticks             = [0.86,0.93,1];
scatter(Wi_v(end)-(Wi_v(end)-Wi_v(end-1))/2,alp_v(end)-...
    (alp_v(end)-alp_v(end-1))/2,65,'rx');
set(gca,'FontSize',15);
axis square;
box on;
ax = gca;
ax.Layer = 'top';
xlabel('Wi','Interpreter','latex');
ylabel('$\alpha$','Interpreter','latex');
clim([0.86,1])
xlim([Wi_v(2),Wi_v(end)])
ylim([alp_v(2),alp_v(end)])
xticks([0,1,2,3,4,5])
yticks([0,0.1,0.2,0.3,0.4,0.5])

f1.Position = [363,410,855,395];