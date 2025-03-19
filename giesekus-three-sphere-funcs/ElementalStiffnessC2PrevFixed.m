function [I1,I2,I3,I4] = ElementalStiffnessC2PrevFixed(shape,tN11,tN12,...
    tN21,tN22,detJ,w2D,TptimeQuad,FixedTime)

% extract shape functions, weights and Jacobian at gaussian quadrature 
% points for non-stokeslet integrals
shapeGauss = FixedTime.shape;
w2DGauss   = FixedTime.w2D;
detJGauss  = FixedTime.detJ{1}; % Jacobians are the same for each element

% extra tensor components
[Tptime11,Tptime12,Tptime21,Tptime22] = ...
    ExtractTensorComponents(TptimeQuad);

% calculate integrals
I1 = shape.t(:,:)*(tN11.*(detJ.*w2D)) + 1/3*shapeGauss.t(:,:)*...
    (Tptime11.*(detJGauss.*w2DGauss));

I2 = shape.t(:,:)*(tN12.*(detJ.*w2D)) + 1/3*shapeGauss.t(:,:)*...
    (Tptime12.*(detJGauss.*w2DGauss));

I3 = shape.t(:,:)*(tN21.*(detJ.*w2D)) + 1/3*shapeGauss.t(:,:)*...
    (Tptime21.*(detJGauss.*w2DGauss));

I4 = shape.t(:,:)*(tN22.*(detJ.*w2D)) + 1/3*shapeGauss.t(:,:)*...
    (Tptime22.*(detJGauss.*w2DGauss));


end