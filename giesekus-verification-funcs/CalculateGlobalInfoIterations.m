function [w2D,shape,detJ,gdshape] = CalculateGlobalInfoIterations(mesh,...
    nquad)

% setup quadrature rule
[x2D,w2D] = GetGaussianQuad2D(nquad);

% get vector of quad points on the reference element
S = [x2D(:,1); x2D(:,2)];

% get reference element shape functions at quad points
shape = CalculateRefShape(S);

% get reference element velocity shape function derivatives at quad points
dshape = CalculateRefShapeDeriv(S);

% calulate jacobian of transformation to the reference element
[Jref,detJ] = CalculateRefJacobian(mesh);

% calculate global shape function derivatives at quad points 
gdshape = CalculateGlobalVelShapeDeriv(mesh,Jref,dshape);

end