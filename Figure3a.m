% Script to plot Figure 3(a) in the manuscript

%% setup

clear all; close all; clc; warning off;

% parameters
alpha  = 0.1;
G      = 1;

% define function (see eq (3.14) in manuscript)
func = @(t0,Wi) (2*alpha-1)./(2*alpha*Wi^2*G)*log((1-alpha*Wi^2*t0.^2)...
    ./(1-alpha*Wi^2*(t0+G).^2)) - ((alpha-1)*(G+2*t0))./((1-alpha*Wi^2*t0.^2)...
    .*(1-alpha*Wi^2*(t0+G).^2))-1;

%% plot

figure;

% Wi = 1 case
subplot(2,2,1);
Wi    = 1;
xeval = linspace(-12,12,1000000);
yeval = func(xeval,Wi);

for ii=1:length(yeval)
    ind(ii) = isreal(yeval(ii));
end

swt    = diff(ind) ~= 0;
nswt   = sum(swt);
swtind = [find(swt == 1), length(xeval)];

title(['Wi$\,=',num2str(Wi),'$'],'Interpreter','latex')
hold on;
line([min(xeval),max(xeval)],[0,0],'LineStyle','-','LineWidth',1.5,'Color','k')
line([1-G/2,1-G/2],[-2,2],'LineStyle',':','LineWidth',1.5,'Color','k')
ylim([-2,2])
xlim([-12,12])
if length(swtind) == 1
    plot(xeval,yeval,'Color','b','LineWidth',1.5)
else
    for ii=1:nswt
        typ(ii) = ind(swtind(ii));
    end

    if ind(end) == 0
        typ(nswt+1) = 0;
    else
        typ(nswt+1) = 1;
    end

    n = 1;
    for ii=1:nswt+1
        if typ(ii) == 0
            clr = 'r';
        else
            clr = 'b';
        end

        plot(xeval(n:swtind(ii)),yeval(n:swtind(ii)),'Color',clr,'LineWidth',1.5)
        n = swtind(ii)+1;
    end
end
box on;
axis square;
set(gca,'FontSize',14)
xlabel('$\tau_0$','Interpreter','latex');
ylabel('$f(\tau_0)$','Interpreter','latex');

% Wi = 2 case
subplot(2,2,2);
Wi    = 2;
xeval = linspace(-9,9,1000000);
yeval = func(xeval,Wi);

for ii=1:length(yeval)
    ind(ii) = isreal(yeval(ii));
end

swt    = diff(ind) ~= 0;
nswt   = sum(swt);
swtind = [find(swt == 1), length(xeval)];

title(['Wi$\,=',num2str(Wi),'$'],'Interpreter','latex')
hold on;
line([min(xeval),max(xeval)],[0,0],'LineStyle','-','LineWidth',1.5,'Color','k')
line([1-G/2,1-G/2],[-2,2],'LineStyle',':','LineWidth',1.5,'Color','k')
ylim([-2,2])
xlim([-9,9])
if length(swtind) == 1
    plot(xeval,yeval,'Color','b','LineWidth',1.5)
else
    for ii=1:nswt
        typ(ii) = ind(swtind(ii));
    end

    if ind(end) == 0
        typ(nswt+1) = 0;
    else
        typ(nswt+1) = 1;
    end

    n = 1;
    for ii=1:nswt+1
        if typ(ii) == 0
            clr = 'r';
        else
            clr = 'b';
        end

        plot(xeval(n:swtind(ii)),yeval(n:swtind(ii)),'Color',clr,'LineWidth',1.5)
        n = swtind(ii)+1;
    end
end
box on;
axis square;
set(gca,'FontSize',14)
xlabel('$\tau_0$','Interpreter','latex');
ylabel('$f(\tau_0)$','Interpreter','latex');


% Wi = 3 case
subplot(2,2,3);
Wi    = 3;
xeval = linspace(-6,6,1000000);
yeval = func(xeval,Wi);

for ii=1:length(yeval)
    ind(ii) = isreal(yeval(ii));
end

swt    = diff(ind) ~= 0;
nswt   = sum(swt);
swtind = [find(swt == 1), length(xeval)];

title(['Wi$\,=',num2str(Wi),'$'],'Interpreter','latex')
hold on;
line([min(xeval),max(xeval)],[0,0],'LineStyle','-','LineWidth',1.5,'Color','k')
line([1-G/2,1-G/2],[-2,2],'LineStyle',':','LineWidth',1.5,'Color','k')
ylim([-2,2])
xlim([-6,6])
if length(swtind) == 1
    plot(xeval,yeval,'Color','b','LineWidth',1.5)
else
    for ii=1:nswt
        typ(ii) = ind(swtind(ii));
    end

    if ind(end) == 0
        typ(nswt+1) = 0;
    else
        typ(nswt+1) = 1;
    end

    n = 1;
    for ii=1:nswt+1
        if typ(ii) == 0
            clr = 'r';
        else
            clr = 'b';
        end

        plot(xeval(n:swtind(ii)),yeval(n:swtind(ii)),'Color',clr,'LineWidth',1.5)
        n = swtind(ii)+1;
    end
end
box on;
axis square;
set(gca,'FontSize',14)
xlabel('$\tau_0$','Interpreter','latex');
ylabel('$f(\tau_0)$','Interpreter','latex');


% Wi = 5 case
subplot(2,2,4);
Wi    = 5;
xeval = linspace(-3,3,1000000);
yeval = func(xeval,Wi);

for ii=1:length(yeval)
    ind(ii) = isreal(yeval(ii));
end

swt    = diff(ind) ~= 0;
nswt   = sum(swt);
swtind = [find(swt == 1), length(xeval)];

title(['Wi$\,=',num2str(Wi),'$'],'Interpreter','latex')
hold on;
line([min(xeval),max(xeval)],[0,0],'LineStyle','-','LineWidth',1.5,'Color','k')
line([1-G/2,1-G/2],[-2,2],'LineStyle',':','LineWidth',1.5,'Color','k')
ylim([-2,2])
xlim([-3,3])
if length(swtind) == 1
    plot(xeval,yeval,'Color','b','LineWidth',1.5)
else
    for ii=1:nswt
        typ(ii) = ind(swtind(ii));
    end

    if ind(end) == 0
        typ(nswt+1) = 0;
    else
        typ(nswt+1) = 1;
    end

    n = 1;
    for ii=1:nswt+1
        if typ(ii) == 0
            clr = 'r';
        else
            clr = 'b';
        end

        plot(xeval(n:swtind(ii)),yeval(n:swtind(ii)),'Color',clr,'LineWidth',1.5)
        n = swtind(ii)+1;
    end
end
box on;
axis square;
set(gca,'FontSize',14)
xlabel('$\tau_0$','Interpreter','latex');
ylabel('$f(\tau_0)$','Interpreter','latex');





