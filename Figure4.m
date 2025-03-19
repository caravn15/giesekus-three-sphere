% Script to plot Figure 4 in the manuscript
% Will take around 10 mins to run

clear all; close all; clc; 
[p,n,e] = fileparts(mfilename('fullpath'));
addpath([p filesep 'giesekus-verification-funcs'],'-end');
rmpath([p filesep 'giesekus-three-sphere-funcs'],'-end');

%% parameters

% domain 
dom.lx = [0,10];
dom.ly = [0,1];

% numerical
dom.nx = 14;
dom.ny = 14;
nquad  = 10;
tol    = 1e-4;

% fluid 
params.pin    = 10;
params.pout   = 0;
params.De = 1;

%% run solve

params.alpha  = 0; params.De = 0;
[u11,~,tau11,~,uAn11] = SolveGiesekusFEM(dom,params,nquad,tol);

params.De = [0.5,1];
params.alpha  = 0.1;
[u22,~,tau22,~,uAn22] = SolveGiesekusFEM(dom,params,nquad,tol);

params.De = [0.5,1];
params.alpha  = 0.2;
[u33,~,tau33,~,uAn33] = SolveGiesekusFEM(dom,params,nquad,tol);

params.De = [0.5,0.75,1];
params.alpha  = 0.3;
[u44,pp44,tau44,~,uAn44] = SolveGiesekusFEM(dom,params,nquad,tol);


params.alpha  = 0; params.De = 0;
[u1,pp1,tau1,~,uAn1] = SolveGiesekusFEM(dom,params,nquad,tol);

params.De = 0.4;
params.alpha  = 0.1;
[u2,pp2,tau2,~,uAn2] = SolveGiesekusFEM(dom,params,nquad,tol,u1,pp1,tau1,uAn1);

params.De = 0.8;
params.alpha  = 0.1;
[u3,pp3,tau3,~,uAn3] = SolveGiesekusFEM(dom,params,nquad,tol,u2,pp2,tau2,uAn2);

params.De = [0.8,1.2];
params.alpha  = 0.1;
[u4,~,tau4,mesh,uAn4] = SolveGiesekusFEM(dom,params,nquad,tol,u3,pp3,tau3,uAn3);


%% plot 

figure;
subplot(1,2,1)
hold on;
plot(uAn11,linspace(0,1,length(uAn11)),'LineWidth',1.5,'Color','b');
plot(uAn22,linspace(0,1,length(uAn11)),'LineWidth',1.5,'Color','b');
plot(uAn33,linspace(0,1,length(uAn11)),'LineWidth',1.5,'Color','b');
plot(uAn44,linspace(0,1,length(uAn11)),'LineWidth',1.5,'Color','b');

y = linspace(0,1,40);
x = 5*ones(size(y));

U      = GetFemVelocityAtPoint(mesh,u11,[x';y']);
[u11,~] = ExtractVectorComponents(U);
plot(u11,y,'o','MarkerEdgeColor','r','LineWidth',1.5,'MarkerFaceColor','r','MarkerSize',4);

U      = GetFemVelocityAtPoint(mesh,u22,[x';y']);
[u22,~] = ExtractVectorComponents(U);
plot(u22,y,'o','MarkerEdgeColor','r','LineWidth',1.5,'MarkerFaceColor','r','MarkerSize',4);

U      = GetFemVelocityAtPoint(mesh,u33,[x';y']);
[u33,~] = ExtractVectorComponents(U);
plot(u33,y,'o','MarkerEdgeColor','r','LineWidth',1.5,'MarkerFaceColor','r','MarkerSize',4);

U      = GetFemVelocityAtPoint(mesh,u44,[x';y']);
[u44,~] = ExtractVectorComponents(U);
plot(u44,y,'o','MarkerEdgeColor','r','LineWidth',1.5,'MarkerFaceColor','r','MarkerSize',4);

box on;
set(gca,'FontSize',13)
xlabel('$u(y)$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
xlim([0,1])


subplot(1,2,2)
hold on;
plot(uAn1,linspace(0,1,length(uAn1)),'LineWidth',1.5,'Color','b');
plot(uAn2,linspace(0,1,length(uAn1)),'LineWidth',1.5,'Color','b');
plot(uAn3,linspace(0,1,length(uAn1)),'LineWidth',1.5,'Color','b');
plot(uAn4,linspace(0,1,length(uAn1)),'LineWidth',1.5,'Color','b');

y = linspace(0,1,40);
x = 5*ones(size(y));

U      = GetFemVelocityAtPoint(mesh,u1,[x';y']);
[u1,~] = ExtractVectorComponents(U);
plot(u1,y,'o','MarkerEdgeColor','r','LineWidth',1.5,'MarkerFaceColor','r','MarkerSize',4);

U      = GetFemVelocityAtPoint(mesh,u2,[x';y']);
[u2,~] = ExtractVectorComponents(U);
plot(u2,y,'o','MarkerEdgeColor','r','LineWidth',1.5,'MarkerFaceColor','r','MarkerSize',4);

U      = GetFemVelocityAtPoint(mesh,u3,[x';y']);
[u3,~] = ExtractVectorComponents(U);
plot(u3,y,'o','MarkerEdgeColor','r','LineWidth',1.5,'MarkerFaceColor','r','MarkerSize',4);

U      = GetFemVelocityAtPoint(mesh,u4,[x';y']);
[u4,~] = ExtractVectorComponents(U);
plot(u4,y,'o','MarkerEdgeColor','r','LineWidth',1.5,'MarkerFaceColor','r','MarkerSize',4);

box on;
set(gca,'FontSize',13)
xlabel('$u(y)$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
xlim([0,1])
