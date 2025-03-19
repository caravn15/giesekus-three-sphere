function PlotStress(Ttime_ind,X,ft,Ttime,refFactor,cbarlim)

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

% create more refined mesh to calculate stress over
xvf       = linspace(lx(1),lx(2),refFactor*2*nx+1);
yvf       = linspace(ly(1),ly(2),refFactor*2*ny+1);
[xvf,yvf] = meshgrid(xvf,yvf);
xvfv = xvf(:); yvfv = yvf(:);

% find position/force index to match Ttime_ind
X_ind_v = 1:nt;
Ttime_ind_v = 1:25:nt;
X_ind = X_ind_v(Ttime_ind_v(Ttime_ind));

% interpolate stress velocity
[T11,T12,T21,T22] = ExtractTensorComponents(Ttime(:,Ttime_ind));
InterpT11 = scatteredInterpolant(xv(:),yv(:),T11);
InterpT12 = scatteredInterpolant(xv(:),yv(:),T12);
InterpT21 = scatteredInterpolant(xv(:),yv(:),T21);
InterpT22 = scatteredInterpolant(xv(:),yv(:),T22);

for ii=1:length(xvf(:))
    T11x0(ii,1) = InterpT11(xvfv(ii),yvfv(ii));
    T12x0(ii,1) = InterpT12(xvfv(ii),yvfv(ii));
    T21x0(ii,1) = InterpT21(xvfv(ii),yvfv(ii));
    T22x0(ii,1) = InterpT22(xvfv(ii),yvfv(ii));
end

% calculate stokeslet tau
for ii=1:3
    dS    = RegStokesletDerivVelocity([xvf(:);yvf(:)],X(:,ii,X_ind),ep);
    dS    = cell2mat(dS);
    dUNdx(:,ii) = dS(:,:,1)*ft(:,ii,X_ind);
    dUNdy(:,ii) = dS(:,:,2)*ft(:,ii,X_ind);

    [duNdx(:,ii),dvNdx(:,ii)] = ExtractVectorComponents(dUNdx(:,ii));
    [duNdy(:,ii),dvNdy(:,ii)] = ExtractVectorComponents(dUNdy(:,ii));

    tN11(:,ii) =  (duNdx(:,ii)+duNdx(:,ii));
    tN12(:,ii) = (duNdy(:,ii)+dvNdx(:,ii));
    tN21(:,ii) = (dvNdx(:,ii)+duNdy(:,ii));
    tN22(:,ii) =  (dvNdy(:,ii)+dvNdy(:,ii));

end
tN11 = sum(tN11,2);
tN12 = sum(tN12,2);
tN21 = sum(tN21,2);
tN22 = sum(tN22,2);

t11 = T11x0 + tN11;
t12 = T12x0 + tN12;
t21 = T21x0 + tN21;
t22 = T22x0 + tN22;
t11 = reshape(t11,size(xvf));
t12 = reshape(t12,size(xvf));
t21 = reshape(t21,size(xvf));
t22 = reshape(t22,size(xvf));

magstress(:,:) = (t11.^2+t22.^2+t12.^2+t21.^2).^(1/2);


hold on;
pcol = pcolor(xvf,yvf,magstress(:,:));

set(pcol,'EdgeColor','None'); shading interp;
colormap(hsv)
line([X(1,1,X_ind),X(1,2,X_ind)],[X(2,1,1),X(2,2,X_ind)],'LineWidth',2,'Color',[0,0,0])
line([X(1,2,X_ind),X(1,3,X_ind)],[X(2,2,1),X(2,3,X_ind)],'LineWidth',2,'Color',[0,0,0])
scatter(X(1,1,X_ind),X(2,1,X_ind),'MarkerEdgeColor',[0,0,0],'MarkerFaceColor',[0,0,0],...
    'SizeData',70);
scatter(X(1,2,X_ind),X(2,2,X_ind),'MarkerEdgeColor',[0,0,0],'MarkerFaceColor',[0,0,0],...
    'SizeData',70);
scatter(X(1,3,X_ind),X(2,3,X_ind),'MarkerEdgeColor',[0,0,0],'MarkerFaceColor',[0,0,0],...
    'SizeData',70);
cb =colorbar;
cb.Location = 'southoutside';
if ~isempty(cbarlim)
clim(cbarlim)
end
xlim([X(1,2,X_ind)-0.3,X(1,2,X_ind)+0.3])
ylim([0.9,1.1])
set(gca,'XTick',[])
set(gca,'YTick',[])
box on;
ax = gca;
ax.Layer = 'top';
set(gca,'FontSize',15);
cb.Label.Interpreter = 'latex';
cb.Label.FontSize    = 15;

end