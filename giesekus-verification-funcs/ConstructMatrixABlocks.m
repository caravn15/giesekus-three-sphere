function AMats = ConstructMatrixABlocks(mesh,gdshape,shape,detJ,w2D)

% initilise full matrices
AMats.A1 = sparse(mesh.np,mesh.nv);
AMats.A2 = sparse(mesh.np,mesh.nv);

for e=1:mesh.nelem
    
    % calculate integrals for each element
    [I1,I2] = ElementalStiffnessA(gdshape{e},shape,detJ{e},w2D);
    
    % assemble full matrix
    AMats.A1(mesh.pCon(e,:), mesh.vCon(e,:)) = ...
        AMats.A1(mesh.pCon(e,:), mesh.vCon(e,:)) + I1;
    
    AMats.A2(mesh.pCon(e,:), mesh.vCon(e,:)) = ...
        AMats.A2(mesh.pCon(e,:), mesh.vCon(e,:)) + I2;
    
end

end