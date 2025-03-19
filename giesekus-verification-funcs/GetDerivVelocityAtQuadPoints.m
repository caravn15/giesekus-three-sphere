function dUQuad = GetDerivVelocityAtQuadPoints(mesh,gdshape,upr)

% extract velocity components
[U1,U2] = ExtractVectorComponents(upr);

for e=1:mesh.nelem
    
    % get velocity components on element
    U1e = U1(mesh.vCon(e,:));
    U2e = U2(mesh.vCon(e,:));
    
    % get velocity x derivatives
    dU1dxQuad = U1e'*squeeze(gdshape{e}(1,:,:));
    dU2dxQuad = U2e'*squeeze(gdshape{e}(1,:,:));
    
    % get velocity y derivatives
    dU1dyQuad = U1e'*squeeze(gdshape{e}(2,:,:));
    dU2dyQuad = U2e'*squeeze(gdshape{e}(2,:,:));
    
    % sort output
    dUQuad{e}.x = [dU1dxQuad(:); dU2dxQuad(:)];
    dUQuad{e}.y = [dU1dyQuad(:); dU2dyQuad(:)];
    
end

end