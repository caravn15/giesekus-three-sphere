function [U,P,T] = SolveNewtonianProblem(model,FixedTime,FixedIts)

%% unpack matrix blocks fixed across time

AMats  = FixedTime.AMats;
BMats  = FixedTime.BMats;
C0Mats = FixedTime.C0Mats;
DVecs  = FixedIts.DVecs;
bc     = FixedIts.bc;

%% assemble full matrix

% zero blocks
Zvv = zeros(2*model.mesh.nv,2*model.mesh.nv);
Zpp = zeros(model.mesh.np,model.mesh.np);
Zpt = zeros(model.mesh.np,4*model.mesh.nt);

% A block
A  = [AMats.A1 AMats.A2];

% B block
B = [BMats.B1x11 BMats.B1x12 BMats.B1x21 BMats.B1x22 ;
     BMats.B2x11 BMats.B2x12 BMats.B2x21 BMats.B2x22];

% F block
C0 = [C0Mats.C011x11 C0Mats.C011x12 C0Mats.C011x21 C0Mats.C011x22 ;
      C0Mats.C012x11 C0Mats.C012x12 C0Mats.C012x21 C0Mats.C012x22 ;
      C0Mats.C021x11 C0Mats.C021x12 C0Mats.C021x21 C0Mats.C021x22 ;
      C0Mats.C022x11 C0Mats.C022x12 C0Mats.C022x21 C0Mats.C022x22];

% assemble matrix
J = [Zvv  -A'   B  ;
     A    Zpp  Zpt;
    -2*B' Zpt' C0 ];

%% assemble full right-hand-side vector

% D vector
D = [DVecs.D1; DVecs.D2];

% assemble vector
a = zeros(2*model.mesh.nv+model.mesh.np+4*model.mesh.nt,1);
a(1:2*model.mesh.nv) = -D;

%% boundary conditions

% apply boundary conditions to matrix
J = ApplyBCsJacobian(model.mesh,J);

% apply boundary conditions to rhs vector
a = ApplyBCsResidual(model.mesh,a,bc);

%% fix pressure

Jz = zeros(1,size(J,2));
Jz(2*model.mesh.nv+1) = 1;
J = [J; Jz];
a = [a;0];

%% solve

% solve system
d = J\a;

% extract solution vectors
U = d(1:2*model.mesh.nv);
P = d(2*model.mesh.nv+1:2*model.mesh.nv+model.mesh.np);
T = d(2*model.mesh.nv+model.mesh.np+1:end);

end