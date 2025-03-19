function PlotVelocity(Utime_ind,X,ft,Utime,refFactor,plotStreams,Ulim)

% unpack parameters
ep = 0.01;
nt = 2501;
lx = [0,2];
ly = [0,2];
nx = 14;
ny = 14;

% recreate mesh
xv      = linspace(lx(1),lx(2),2*nx+1);
yv      = linspace(ly(1),ly(2),2*ny+1);
[xv,yv] = meshgrid(xv,yv);

% create more refined mesh to calculate velocities over
xvf       = linspace(lx(1),lx(2),refFactor*2*nx+1);
yvf       = linspace(ly(1),ly(2),refFactor*2*ny+1);
[xvf,yvf] = meshgrid(xvf,yvf);
xvfv = xvf(:); yvfv = yvf(:);

% find position/force index to match Utime_ind
X_ind_v = 1:nt;
Utime_ind_v = 1:25:nt;
X_ind = X_ind_v(Utime_ind_v(Utime_ind));

% interpolate fem velocity
[U,V] = ExtractVectorComponents(Utime(:,Utime_ind));
InterpU = scatteredInterpolant(xv(:),yv(:),U);
InterpV = scatteredInterpolant(xv(:),yv(:),V);
for ii=1:length(xvf(:))
    Ux0(ii,1) = InterpU(xvfv(ii),yvfv(ii));
    Vx0(ii,1) = InterpV(xvfv(ii),yvfv(ii));
end

% calculate stokeslet velocity
for ii=1:3
    S = RegStokesletVelocity([xvf(:);yvf(:)],X(:,ii,X_ind),ep);
    UN = S*ft(:,ii,X_ind);
    [uN(:,ii),vN(:,ii)] = ExtractVectorComponents(UN);
end
uN = sum(uN,2);
vN = sum(vN,2);

% calculate full velocity for plot
u = uN + Ux0;
v = vN + Vx0;
ux = reshape(u,size(xvf));
uy = reshape(v,size(yvf));
Uflow = ((ux.^2 + uy.^2).^(0.5));

% plot
hold on;
pcol = pcolor(xvf,yvf,Uflow);
set(pcol,'EdgeColor','None');
shading interp;
cb = colorbar;
cb.Label.Interpreter = 'latex';
cb.Label.FontSize    = 15;
set(gca,'FontSize',15);
axis square;
box on;
ax = gca;
ax.Layer = 'top';

if plotStreams == 1
    stream = streamslice(xvf,yvf,ux,uy,0.6);
    set(stream,'color',[163,198,255]./255,'LineWidth',0.8)
end

line([X(1,1,X_ind),X(1,2,X_ind)],[X(2,1,1),X(2,2,X_ind)],'LineWidth',2,'Color','red')
line([X(1,2,X_ind),X(1,3,X_ind)],[X(2,2,1),X(2,3,X_ind)],'LineWidth',2,'Color','red')
scatter(X(1,1,X_ind),X(2,1,X_ind),'MarkerEdgeColor','r','MarkerFaceColor','r',...
    'SizeData',70);
scatter(X(1,2,X_ind),X(2,2,X_ind),'MarkerEdgeColor','r','MarkerFaceColor','r',...
    'SizeData',70);
scatter(X(1,3,X_ind),X(2,3,X_ind),'MarkerEdgeColor','r','MarkerFaceColor','r',...
    'SizeData',70);

if ~isempty(Ulim)
    clim(Ulim);
end
clim([0,0.11])
colorbar off;
set(gca,'XTick',[])
set(gca,'YTick',[])

box on;
ax = gca;
ax.Layer = 'top';
end