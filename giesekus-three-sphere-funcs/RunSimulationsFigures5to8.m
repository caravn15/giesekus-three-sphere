function RunSimulationsFigures5to8()

% this will take a (very) long time and I don't recommend running it on 
% a single computer!

close all; clear all; clc;

% parameters
chi    = pi/2;
Amp    = 0.8;
ep     = 0.01;
dt     = 0.01;
tmin   = 0;
kappa  = 0.8;
nbeats = 25;
lr     = 0.2;
ng     = 14;

% swimmer initial condition
X0(:,1) = [1-lr;1];
X0(:,2) = [1;1];
X0(:,3) = [1+lr;1];

% define alph and Wi
alp_v        = linspace(0,0.5,25);
Wi_v         = linspace(0,5,25);
[Wi_M,alp_M] = meshgrid(Wi_v,alp_v);
alp_V        = alp_M(:);
Wi_V         = Wi_M(:);

for jj=1:5 % solves problem in 5 sets of 25 beats 

    if jj~=1

        global FixedItsPrev;

        tmin            = model.tmax;
        X0              = Xt(:,:,end-1);
        dataprev.Tptime = Tptime;
        dataprev.Uptime = Uptime;
        dataprev.Pptime = Pptime;
        FixedItsPrev    = FinalFixedIts;

    end

    for kk = 1:length(alp_V)

        alph = alp_V(kk);
        Wi   = Wi_V(kk);

        if jj  == 1
            [model,Xt,ft,velX,Utime,Ptime,Ttime,Uptime,Pptime,Tptime,FinalFixedIts,vels] = ...
                RunGiesekusThreeSphereSwimmer(X0,kappa,Amp,chi,lr,ep,alph,...
                Wi,dt,nbeats,tmin,ng);
        else
            [model,Xt,ft,velX,Utime,Ptime,Ttime,Uptime,Pptime,Tptime,FinalFixedIts,vels] = ...
                RunGiesekusThreeSphereSwimmer(X0,kappa,Amp,chi,lr,ep,alph,...
                Wi,dt,nbeats,tmin,ng,dataprev);
        end

        clear dataprev FixedItsPrev;

        save(['data/set_',num2str(jj),'_sim_',num2str(kk-1),'.mat']);
    end
end

end