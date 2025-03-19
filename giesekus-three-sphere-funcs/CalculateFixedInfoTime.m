function FixedTime = CalculateFixedInfoTime(model)

%% setup basis functions and integral mappings

% setup quadrature rules
[x2D,w2D] = GetGaussianQuad2D(model.num.gaussquad);

% get vector of quad points on the reference element
S = [x2D(:,1); x2D(:,2)];

% get reference element shape functions at quad points
[shape] = CalculateRefShape(S);

% get reference element velocity shape function derivatives at quad points
dshape = CalculateRefShapeDeriv(S);

% calulate jacobian of transformation to the reference element
[Jref,detJ] = CalculateRefJacobian(model.mesh);

% calculate global shape function derivatives at quad points 
gdshape = CalculateGlobalVelShapeDeriv(model.mesh,Jref,dshape);

%% construct matrix blocks

% momentum equation (BMats gives T^(m,r):D(v), P^(m,r)(nabla.v) is 
% transpose of AMats)
BMats = ConstructMatrixBBlocks(model.mesh,gdshape,shape,detJ,w2D); 

% continuity equation (AMats gives q(nabla.U^(m,r)))
AMats = ConstructMatrixABlocks(model.mesh,gdshape,shape,detJ,w2D);

% stress constitutive equation (C0Mats gives T^(m,r):theta, 
% D(U^(m,r)):theta is transpose of BMats)
C0Mats = ConstructMatrixC0Blocks(model.mesh,shape,detJ,w2D);

%% package into structure

FixedTime.BMats   = BMats;
FixedTime.AMats   = AMats;
FixedTime.C0Mats  = C0Mats;
FixedTime.shape   = shape;
FixedTime.gdshape = gdshape;
FixedTime.detJ    = detJ;
FixedTime.w2D     = w2D;

end