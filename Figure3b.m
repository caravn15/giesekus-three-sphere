% Script to plot Figure 3(b) in the manuscript

%% setup

clear all; close all; clc; warning off;

% parameters
alpha  = 0.1;
G      = 1;

%% plot

figure;

Wi    = 2;
func  = @(t0) (2*alpha-1)./(2*alpha*Wi^2*G)*log((1-alpha*Wi^2*t0.^2)...
    ./(1-alpha*Wi^2*(t0+G).^2)) - ((alpha-1)*(G+2*t0))./((1-alpha*Wi^2*t0.^2)...
    .*(1-alpha*Wi^2*(t0+G).^2))-1;
t0_lower = fsolve(func,-0.5);
t0_upper = fsolve(func,3);

vel = @(y) -1/(2*alpha*G*Wi^2)*((2*(alpha-1))./...
    (1-alpha*Wi^2*(t0_upper+G*y).^2)+(2*alpha-1)*...
    log(1-alpha*Wi^2*(t0_upper+G*y).^2));
C   = -vel(0);
y   = linspace(0,1,100);
u_upper = vel(y) + C;

subplot(2,2,1);
hold on;
% line([0,0],[0,1],'LineWidth',1,'Color','k')
plot(u_upper,y,'k','LineWidth',1.5)
set(gca,'FontSize',14)
box on;
xlabel('$u(y)$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
xlim([-0.4,1])

vel = @(y) -1/(2*alpha*G*Wi^2)*((2*(alpha-1))./...
    (1-alpha*Wi^2*(t0_lower+G*y).^2)+(2*alpha-1)*...
    log(1-alpha*Wi^2*(t0_lower+G*y).^2));
C   = -vel(0);
y   = linspace(0,1,100);
u_lower = vel(y) + C;

subplot(2,2,2);
hold on;
% line([0,0],[0,1],'LineWidth',1,'Color','k')
plot(u_lower,y,'k','LineWidth',1.5)
set(gca,'FontSize',14)
box on;
xlabel('$u(y)$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
xlim([-0.4,1])

Wi    = 5;
func  = @(t0) (2*alpha-1)./(2*alpha*Wi^2*G)*log((1-alpha*Wi^2*t0.^2)...
    ./(1-alpha*Wi^2*(t0+G).^2)) - ((alpha-1)*(G+2*t0))./((1-alpha*Wi^2*t0.^2)...
    .*(1-alpha*Wi^2*(t0+G).^2))-1;
t0_lower = fsolve(func,-0.5);
t0_upper = fsolve(func,1);

vel = @(y) -1/(2*alpha*G*Wi^2)*((2*(alpha-1))./...
    (1-alpha*Wi^2*(t0_upper+G*y).^2)+(2*alpha-1)*...
    log(1-alpha*Wi^2*(t0_upper+G*y).^2));
C   = -vel(0);
y   = linspace(0,1,100);
u_upper = vel(y) + C;

subplot(2,2,3);
hold on;
% line([0,0],[0,1],'LineWidth',1,'Color','k')
plot(u_upper,y,'k','LineWidth',1.5)
set(gca,'FontSize',14)
box on;
xlabel('$u(y)$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
xlim([-0.4,1])

vel = @(y) -1/(2*alpha*G*Wi^2)*((2*(alpha-1))./...
    (1-alpha*Wi^2*(t0_lower+G*y).^2)+(2*alpha-1)*...
    log(1-alpha*Wi^2*(t0_lower+G*y).^2));
C   = -vel(0);
y   = linspace(0,1,100);
u_lower = vel(y) + C;

subplot(2,2,4);
hold on;
% line([0,0],[0,1],'LineWidth',1,'Color','k')
plot(u_lower,y,'k','LineWidth',1.5)
set(gca,'FontSize',14)
box on;
xlabel('$u(y)$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
xlim([-0.4,1])