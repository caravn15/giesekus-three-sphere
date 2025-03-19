function [u0,p0,tau0,uAn] = SolveNewtonianProblem(mesh,params,w2D,shape,...
    gdshape,detJ)

%% calculate analytic solution

uAn = CalculateAnalyticSolution(params,mesh.lx(2),...
    mesh.yv(mesh.vbind{1}),'Newtonian');

%% construct matrix blocks

% momentum equation
BMats = ConstructMatrixBBlocks(mesh,gdshape,shape,detJ,w2D);

% continuity equation 
AMats = ConstructMatrixABlocks(mesh,gdshape,shape,detJ,w2D);

% stress constitutive equation
FMats = ConstructMatrixFBlocks(mesh,shape,detJ,w2D);

%% assemble full matrix

% zero blocks
Zvv = zeros(2*mesh.nv,2*mesh.nv);
Zpp = zeros(mesh.np,mesh.np);
Zpt = zeros(mesh.np,4*mesh.nt);

% A block
A  = [AMats.A1 AMats.A2];

% B block
B = [BMats.B1x11 BMats.B1x12 BMats.B1x21 BMats.B1x22 ;
     BMats.B2x11 BMats.B2x12 BMats.B2x21 BMats.B2x22];

% F block
F = [FMats.F11x11 FMats.F11x12 FMats.F11x21 FMats.F11x22 ;
     FMats.F12x11 FMats.F12x12 FMats.F12x21 FMats.F12x22 ;
     FMats.F21x11 FMats.F21x12 FMats.F21x21 FMats.F21x22 ;
     FMats.F22x11 FMats.F22x12 FMats.F22x21 FMats.F22x22];

% assemble matrix
J = [Zvv  -A'   B  ;
     A    Zpp  Zpt;
    -2*B' Zpt' F ];

%% boundary conditions

% apply boundary conditions to matrix
J = ApplyBCsJacobian(mesh,J);

% apply boundary conditions to rhs vector
a = zeros(2*mesh.nv+mesh.np+4*mesh.nt,1);
a = ApplyBCsResidual(mesh,a,uAn);

%% solve

% solve system
d = J\a;

% extract solution vectors
u0   = d(1:2*mesh.nv);
p0   = d(2*mesh.nv+1:2*mesh.nv+mesh.np);
tau0 = d(2*mesh.nv+mesh.np+1:end);

end