function mesh = GenerateMeshRect(model)

% unpack structures
lx = model.mesh.lx; ly = model.mesh.ly; 
nx = model.mesh.nx; ny = model.mesh.ny; 

% number of elements
mesh.nelem = nx*ny;

% other parameters
mesh.lx = lx; mesh.ly = ly;
mesh.nx = nx; mesh.ny = ny;

% pressure nodes
xp      = linspace(lx(1),lx(2),nx+1);
yp      = linspace(ly(1),ly(2),ny+1);
[xp,yp] = meshgrid(xp,yp);
mesh.xp = xp(:);
mesh.yp = yp(:);
mesh.np = length(mesh.xp);

% pressure boundary nodes
mesh.pbind{1} = 1:ny+1;
mesh.pbind{2} = ny+1:ny+1:mesh.np;
mesh.pbind{3} = mesh.np-ny:mesh.np;
mesh.pbind{4} = 1:ny+1:mesh.np-ny;

% velocity and stress nodes
xv      = linspace(lx(1),lx(2),2*nx+1);
yv      = linspace(ly(1),ly(2),2*ny+1);
[xv,yv] = meshgrid(xv,yv);
mesh.xv = xv(:);
mesh.yv = yv(:);
mesh.nv = length(mesh.xv);
mesh.xt = xv(:);
mesh.yt = yv(:);
mesh.nt = length(mesh.xt);

% velocity and stress boundary nodes
mesh.vbind{1} = 1:2*ny+1;
mesh.vbind{2} = 2*ny+1:2*ny+1:mesh.nv;
mesh.vbind{3} = mesh.nv-2*ny:mesh.nv;
mesh.vbind{4} = 1:2*ny+1:mesh.nv-2*ny;
mesh.tbind{1} = 1:2*ny+1;
mesh.tbind{2} = 2*ny+1:2*ny+1:mesh.nv;
mesh.tbind{3} = mesh.nv-2*ny:mesh.nv;
mesh.tbind{4} = 1:2*ny+1:mesh.nv-2*ny;

% pressure connectivity
count = 1; pCon = zeros(mesh.nelem,4);
for jj=1:nx
    for ii=1:ny
        pCon(count,:) = [ii+(jj-1)*(ny+1), ii+(ny+1)+(jj-1)*(ny+1), ...
            ii+1+(ny+1)+(jj-1)*(ny+1), ii+1+(jj-1)*(ny+1)];
        count = count + 1;
    end
end
mesh.pCon = pCon;

% velocity and stress connectivity
count = 1; vCon = zeros(mesh.nelem,9);
for jj=1:nx
    for ii=1:2:2*ny
        vCon(count,:) = [ii+(jj-1)*2*(2*ny+1), ...
            ii+(jj-1)*2*(2*ny+1)+(2*ny+1), ...
            ii+(jj-1)*2*(2*ny+1)+2*(2*ny+1), ...
            ii+1+(jj-1)*2*(2*ny+1)+2*(2*ny+1), ...
            ii+2+(jj-1)*2*(2*ny+1)+2*(2*ny+1), ...
            ii+2+(jj-1)*2*(2*ny+1)+(2*ny+1), ...
            ii+2+(jj-1)*2*(2*ny+1), ii+1+(jj-1)*2*(2*ny+1), ...
            ii+1+(jj-1)*2*(2*ny+1)+(2*ny+1)];
        count = count + 1;
    end
end
mesh.vCon = vCon;
mesh.tCon = vCon;
mesh.elemout = mesh.nelem-ny+1:mesh.nelem;
mesh.elemin = 1:ny;

% calculate limits of each element 
mesh.xlim = zeros(mesh.nelem,2);
mesh.ylim = zeros(mesh.nelem,2);

for e=1:mesh.nelem
    
    xe = mesh.xv(mesh.vCon(e,:));
    ye = mesh.yv(mesh.vCon(e,:));
    
    mesh.xlim(e,:) = [min(xe),max(xe)];
    mesh.ylim(e,:) = [min(ye),max(ye)];

end

end