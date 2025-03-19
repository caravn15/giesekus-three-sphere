function RunSimulationsFigure2a_Newt(npdir)

close all; clc;

% parameters
chi    = pi/2;
Amp    = 0.8;
ep     = 0.01;
Wi     = 0;
alph   = 0;
dt     = 0.01;
tmin   = 0;
kappa  = 0.8;
nbeats = 0;
lr     = 0.2;

% swimmer initial condition
X0(:,1) = [1-lr;1];
X0(:,2) = [1;1];
X0(:,3) = [1+lr;1];


%% run simulations

for ii=2:npdir

    [~,Xt,ft,~,~,~,~,~,~,~,~,vels] = RunGiesekusThreeSphereSwimmer(...
        X0,kappa,Amp,chi,lr,ep,alph,Wi,dt,nbeats,tmin,ii);

    % save data
    save(['data/data_fig_2a_newt_',num2str(ii),'.mat']);

end

end