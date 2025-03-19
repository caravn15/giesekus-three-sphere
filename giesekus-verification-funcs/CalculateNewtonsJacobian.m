function J = CalculateNewtonsJacobian(mesh,params,dUQuad,TQuad,w2D,...
    shape,gdshape,detJ)

%% construct matrix blocks

% momentum equation
BMats = ConstructMatrixBBlocks(mesh,gdshape,shape,detJ,w2D);

% continuity equation 
AMats = ConstructMatrixABlocks(mesh,gdshape,shape,detJ,w2D);

% stress constitutive equation
CMats = ConstructMatrixCBlocks(mesh,gdshape,shape,TQuad,detJ,w2D);
DMats = ConstructMatrixDBlocks(mesh,shape,TQuad,detJ,w2D);
% EMats = ConstructMatrixEBlocks(mesh,shape,dUQuad,detJ,w2D);
FMats = ConstructMatrixFBlocks(mesh,shape,detJ,w2D);

% construct full matrix
J = AssembleJacobianBlocks(mesh,params,AMats,BMats,CMats,DMats,FMats);

% apply boundary conditions
J = ApplyBCsJacobian(mesh,J);

end