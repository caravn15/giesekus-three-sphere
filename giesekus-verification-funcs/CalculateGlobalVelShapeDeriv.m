function gdshape = CalculateGlobalVelShapeDeriv(mesh,Jref,dshape)

gdshape = cell(mesh.nelem,1);

for e=1:mesh.nelem
    
    gdshape{e} = zeros(size(dshape.v));
    
    % calculate inverse of jacobian
    invJ = inv(Jref{e});
    
    % calculate global shape function derivatives
    gdshape{e}(1,:,:) = invJ(1,1)*dshape.v(1,:,:);
    gdshape{e}(2,:,:) = invJ(2,2)*dshape.v(2,:,:);
end

end