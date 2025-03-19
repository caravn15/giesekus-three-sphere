% Script to plot Figure 2(a) in the manuscript

clear all; close all; clc;
[p,~,~] = fileparts(mfilename('fullpath'));
addpath([p '/giesekus-three-sphere-funcs']);
rmpath([p filesep 'giesekus-verification-funcs'],'-end');

%% setup

clear all; close all; clc;

% parameters
ep    = 0.01;
lx    = [0,2];
ly    = [0,2];
npdir = 30;

% get data
try
    load(['data/data_fig_2a_newt_',num2str(npdir),'.mat']);
    load(['data/data_fig_2a_nonnewt_',num2str(npdir),'.mat']);
catch
    RunSimulationsFigure2a_Newt(npdir);
    RunSimulationsFigure2a_NonNewt(npdir);
end

%% calculations for Newtonian case

% load data
for jj = 2:npdir
    load(['data/data_fig_2a_newt_',num2str(jj),'.mat']);
    vel{jj} = vels{1};
end

% create refined mesh to calculate errors over
xvf       = linspace(lx(1),lx(2),2*npdir+1);
yvf       = linspace(ly(1),ly(2),2*npdir+1);
[xvf,yvf] = meshgrid(xvf,yvf);
xvfv = xvf(:); yvfv = yvf(:);

for jj=2:30

    % create mesh with jj elements per direction
    xv      = linspace(lx(1),lx(2),2*jj+1);
    yv      = linspace(ly(1),ly(2),2*jj+1);
    [xv,yv] = meshgrid(xv,yv);

    % interpolate fem velocity to refined mesh points
    [U,V]   = ExtractVectorComponents(vel{jj}.U);
    InterpU = scatteredInterpolant(xv(:),yv(:),U);
    InterpV = scatteredInterpolant(xv(:),yv(:),V);
    for ii=1:length(xvf(:))
        Ux0(ii,1) = InterpU(xvfv(ii),yvfv(ii));
        Vx0(ii,1) = InterpV(xvfv(ii),yvfv(ii));
    end

    % calculate stokeslet velocities over refined mesh
    for ii=1:3
        S  = RegStokesletVelocity([xvf(:);yvf(:)],Xt(:,ii,1),ep);
        UN = S*ft(:,ii,1);
        [uN(:,ii),vN(:,ii)] = ExtractVectorComponents(UN);
    end
    uN = sum(uN,2);
    vN = sum(vN,2);

    % calculate full velocity for plot
    u(:,jj) = uN + Ux0;
    v(:,jj) = vN + Vx0;

end

% calculate MSE
for jj = 2:npdir
    e1           = (u(:,jj)-u(:,end)).^2;
    e2           = sum(e1);
    err_newt(jj) = e2/length(xvfv);
end

%% calculations for non-Newtonian case

% load data
for jj = 2:npdir
    load(['data/data_fig_2a_nonnewt_',num2str(jj),'.mat']);
    vel{jj} = vels{1};
end

% create refined mesh to calculate errors over
xvf       = linspace(lx(1),lx(2),2*npdir+1);
yvf       = linspace(ly(1),ly(2),2*npdir+1);
[xvf,yvf] = meshgrid(xvf,yvf);
xvfv = xvf(:); yvfv = yvf(:);

for jj=2:30

    % create mesh with jj elements per direction
    xv      = linspace(lx(1),lx(2),2*jj+1);
    yv      = linspace(ly(1),ly(2),2*jj+1);
    [xv,yv] = meshgrid(xv,yv);

    % interpolate fem velocity to refined mesh points
    [U,V]   = ExtractVectorComponents(vel{jj}.U);
    InterpU = scatteredInterpolant(xv(:),yv(:),U);
    InterpV = scatteredInterpolant(xv(:),yv(:),V);
    for ii=1:length(xvf(:))
        Ux0(ii,1) = InterpU(xvfv(ii),yvfv(ii));
        Vx0(ii,1) = InterpV(xvfv(ii),yvfv(ii));
    end

    % calculate stokeslet velocities over refined mesh
    for ii=1:3
        S  = RegStokesletVelocity([xvf(:);yvf(:)],Xt(:,ii,1),ep);
        UN = S*ft(:,ii,1);
        [uN(:,ii),vN(:,ii)] = ExtractVectorComponents(UN);
    end
    uN = sum(uN,2);
    vN = sum(vN,2);

    % calculate full velocity for plot
    u(:,jj) = uN + Ux0;
    v(:,jj) = vN + Vx0;

end

% calculate MSE
for jj = 2:npdir
    e1              = (u(:,jj)-u(:,end)).^2;
    e2              = sum(e1);
    err_nonnewt(jj) = e2/length(xvfv);
end

%% plot

figure;
hold on;
plot((2:npdir).^2,err_newt(2:end),'-o','LineWidth',1.5,...
    'MarkerFaceColor',[0 0.4470 0.7410])
plot((2:npdir).^2,err_nonnewt(2:end),'-^','LineWidth',1.5,...
    'MarkerFaceColor',[0.8500 0.3250 0.0980])
set(gca, 'YScale', 'log');
set(gca,'FontSize',14);
xlim([0,900]);
legend({'Wi $ =0$, $\alpha=0$ (Newtonian)','Wi $ =3$, $\alpha=0.5$ (Giesekus)'},...
    'Interpreter','latex');
grid on;
box on;
xlabel('$N_e$','Interpreter','latex');
ylabel('MSE','Interpreter','latex');