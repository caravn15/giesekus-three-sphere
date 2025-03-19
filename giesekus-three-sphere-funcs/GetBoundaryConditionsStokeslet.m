function bc = GetBoundaryConditionsStokeslet(model,Xt,ft)

% unpack required parameters
ep   = model.ep;
mesh = model.mesh;

% get velocity boundary indices from mesh structure
binds = [mesh.vbind{1}'; mesh.vbind{2}'; mesh.vbind{3}'; mesh.vbind{4}'];
binds = unique(binds);

% get mesh boundary coordinates
X = [mesh.xv(binds,1); mesh.yv(binds,1)];

for ii=1:3
    
% calculate stokeslet velocity solution at boundary nodes
S       = RegStokesletVelocity(X,Xt(:,ii),ep);
UN      = S*ft(:,ii);
[uN(:,ii),vN(:,ii)] = ExtractVectorComponents(UN);

end

% sort output
bc.uN   = sum(uN,2);
bc.vN   = sum(vN,2);
bc.inds = binds;

end