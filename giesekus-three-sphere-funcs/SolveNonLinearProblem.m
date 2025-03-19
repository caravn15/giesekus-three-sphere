function [Un,Pn,Tn] = SolveNonLinearProblem(model,TpIt,FixedTime,FixedIts)

%% setup

global FixedItsPrev;

if isempty(FixedItsPrev) && model.tmin == 0
    calc = 1;
else
    calc = 0;
end

%% unpack information

% unpack information fixed across time
AMats   = FixedTime.AMats;
BMats   = FixedTime.BMats;
C0Mats  = FixedTime.C0Mats;
shape   = FixedTime.shape;
gdshape = FixedTime.gdshape;
detJ    = FixedTime.detJ;
w2D     = FixedTime.w2D;

% unpack information fixed across iterations
bc          = FixedIts.bc;
DVecs       = FixedIts.DVecs;
C1FixedVecs = FixedIts.C1FixedVecs;
H1FixedVecs = FixedIts.H1FixedVecs;
C1FixedMats = FixedIts.C1FixedMats;
H1FixedMats = FixedIts.H1FixedMats;
if calc == 0
    C2CurFixedVecs  = FixedIts.C2CurFixedVecs;
    C2PrevFixedVecs = FixedIts.C2PrevFixedVecs;
    H0FixedVecs     = FixedIts.H0FixedVecs;
    C2CurFixedMats  = FixedIts.C2CurFixedMats;
    H0FixedMats     = FixedIts.H0FixedMats;
end

%% construct matrix blocks - stress constitutive equation

% evaluate stress from previous iteration at quadrature points
TpItQuad  = GetStressAtQuadPoints(model.mesh,shape,TpIt);

% evaluate stress derivatives from previous iteration at quadrature points
dTItQuad = GetFemStressDerivatives(TpIt,model.mesh,gdshape);

% H1Mats gives (T^(m,r-1).nabla(U^(m,r))+(nabla(U^(m,r)))^T.T^(m,r-1)):theta
H1Mats = ConstructMatrixH1Blocks(model.mesh,gdshape,shape,TpItQuad,detJ,w2D);

% C1Mats gives (T^(m,r-1).T^(m,r)):theta
C1Mats = ConstructMatrixC1Blocks(model.mesh,shape,TpItQuad,detJ,w2D);

if calc == 0

    % H0Mats gives U^(m,r).nabla(T^(m,r-1)):theta
    H0Mats = ConstructMatrixH0Blocks(model.mesh,shape,dTItQuad,detJ,w2D);

end

% model parameters
alpha = model.alpha;
Wi    = model.Wi;
dt    = model.dt;

% A block
AM  = [AMats.A1 AMats.A2];

% B block
BM = [BMats.B1x11 BMats.B1x12 BMats.B1x21 BMats.B1x22 ;
    BMats.B2x11 BMats.B2x12 BMats.B2x21 BMats.B2x22];

% H1 block
H1M = [H1FixedMats.H111x1 H1FixedMats.H111x2 ;
    H1FixedMats.H112x1 H1FixedMats.H112x2 ;
    H1FixedMats.H121x1 H1FixedMats.H121x2 ;
    H1FixedMats.H122x1 H1FixedMats.H122x2] + ...
    [H1Mats.H111x1      H1Mats.H111x2 ;
    H1Mats.H112x1      H1Mats.H112x2 ;
    H1Mats.H121x1      H1Mats.H121x2 ;
    H1Mats.H122x1      H1Mats.H122x2];

% C0 block
C0M = [C0Mats.C011x11 C0Mats.C011x12 C0Mats.C011x21 C0Mats.C011x22 ;
    C0Mats.C012x11 C0Mats.C012x12 C0Mats.C012x21 C0Mats.C012x22 ;
    C0Mats.C021x11 C0Mats.C021x12 C0Mats.C021x21 C0Mats.C021x22 ;
    C0Mats.C022x11 C0Mats.C022x12 C0Mats.C022x21 C0Mats.C022x22];

% C1 block
C1M = [C1FixedMats.C111x11 C1FixedMats.C111x12 C1FixedMats.C111x21 C1FixedMats.C111x22 ;
    C1FixedMats.C112x11 C1FixedMats.C112x12 C1FixedMats.C112x21 C1FixedMats.C112x22 ;
    C1FixedMats.C121x11 C1FixedMats.C121x12 C1FixedMats.C121x21 C1FixedMats.C121x22 ;
    C1FixedMats.C122x11 C1FixedMats.C122x12 C1FixedMats.C122x21 C1FixedMats.C122x22] + ...
    [C1Mats.C111x11      C1Mats.C111x12      C1Mats.C111x21      C1Mats.C111x22 ;
    C1Mats.C112x11      C1Mats.C112x12      C1Mats.C112x21      C1Mats.C112x22 ;
    C1Mats.C121x11      C1Mats.C121x12      C1Mats.C121x21      C1Mats.C121x22 ;
    C1Mats.C122x11      C1Mats.C122x12      C1Mats.C122x21      C1Mats.C122x22];

if calc == 0

    % C2 block
    C2M = [C2CurFixedMats.C211x11 C2CurFixedMats.C211x12 C2CurFixedMats.C211x21 C2CurFixedMats.C211x22 ;
        C2CurFixedMats.C212x11 C2CurFixedMats.C212x12 C2CurFixedMats.C212x21 C2CurFixedMats.C212x22 ;
        C2CurFixedMats.C221x11 C2CurFixedMats.C221x12 C2CurFixedMats.C221x21 C2CurFixedMats.C221x22 ;
        C2CurFixedMats.C222x11 C2CurFixedMats.C222x12 C2CurFixedMats.C222x21 C2CurFixedMats.C222x22];

    % H0 block
    H0M = [H0FixedMats.H011x1 H0FixedMats.H011x2 ;
        H0FixedMats.H012x1 H0FixedMats.H012x2 ;
        H0FixedMats.H021x1 H0FixedMats.H021x2 ;
        H0FixedMats.H022x1 H0FixedMats.H022x2] + ...
        [H0Mats.H011x1      H0Mats.H011x2 ;
        H0Mats.H012x1      H0Mats.H012x2 ;
        H0Mats.H021x1      H0Mats.H021x2 ;
        H0Mats.H022x1      H0Mats.H022x2];

else
    H0M = 0;
    C2M = 0;
end

%% construct vector components - constitutive equation

% C1Vecs gives (T^(m,r-1).tau_N^m):theta
C1Vecs = ConstructVectorC1(model,FixedIts,TpIt);

% H1Vecs gives (T^(m,r-1).nabla(u_N^m)+(nabla(u_N^m))^T.T^(m,r-1)):theta
H1Vecs = ConstructVectorH1(model,FixedIts,TpIt);

if calc == 0

    % H0Vecs gives (u_N^m.nabla(T^(m,r-1))):theta
    H0Vecs = ConstructVectorH0(model,FixedIts,TpIt);

end

% D vector
DV = [DVecs.D1; DVecs.D2];

% C1 vector
C1V = [C1Vecs.C11; C1Vecs.C12; C1Vecs.C13; C1Vecs.C14] + ...
    [C1FixedVecs.C11; C1FixedVecs.C12; C1FixedVecs.C13; C1FixedVecs.C14];

% H1 vector
H1V = [H1Vecs.H11; H1Vecs.H12; H1Vecs.H13; H1Vecs.H14] + ...
    [H1FixedVecs.H11; H1FixedVecs.H12; H1FixedVecs.H13; H1FixedVecs.H14];

if calc == 0

    % H0 vector
    H0V = [H0Vecs.H01; H0Vecs.H02; H0Vecs.H03; H0Vecs.H04] + ...
        [H0FixedVecs.H01; H0FixedVecs.H02; H0FixedVecs.H03; H0FixedVecs.H04];

    % C2 vector
    C2V = [C2PrevFixedVecs.C21; C2PrevFixedVecs.C22; C2PrevFixedVecs.C23; C2PrevFixedVecs.C24] - ...
        [C2CurFixedVecs.C21; C2CurFixedVecs.C22; C2CurFixedVecs.C23; C2CurFixedVecs.C24];

else
    H0V = 0;
    C2V = 0;
end

%% assemble system

[J,a] = AssembleSystem(alpha,Wi,dt,model,calc,H0M,H1M,BM,AM,C0M,C1M,C2M,...
    C1V,C2V,H0V,H1V,DV);

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
Un = d(1:2*model.mesh.nv);
Pn = d(2*model.mesh.nv+1:2*model.mesh.nv+model.mesh.np);
Tn = d(2*model.mesh.nv+model.mesh.np+1:end);

end