function CMats = ConstructMatrixCBlocks(mesh,gdshape,shape,TQuad,detJ,w2D)

% initilise full matrices
CMats.C1x11 = sparse(mesh.nt,mesh.nv);
CMats.C1x12 = sparse(mesh.nt,mesh.nv);
CMats.C1x21 = sparse(mesh.nt,mesh.nv);
CMats.C1x22 = sparse(mesh.nt,mesh.nv);

CMats.C2x11 = sparse(mesh.nt,mesh.nv);
CMats.C2x12 = sparse(mesh.nt,mesh.nv);
CMats.C2x21 = sparse(mesh.nt,mesh.nv);
CMats.C2x22 = sparse(mesh.nt,mesh.nv);

for e=1:mesh.nelem
    
    % calculate integrals for each element
    [I1x11,I1x12,I1x21,I1x22,I2x11,I2x12,I2x21,I2x22] = ...
        ElementalStiffnessC(TQuad{e},gdshape{e},shape,detJ{e},w2D);
    
    % assemble full matrix
    CMats.C1x11(mesh.tCon(e,:), mesh.vCon(e,:)) = ...
        CMats.C1x11(mesh.tCon(e,:), mesh.vCon(e,:)) + I1x11;
    
    CMats.C1x12(mesh.tCon(e,:), mesh.vCon(e,:)) = ...
        CMats.C1x12(mesh.tCon(e,:), mesh.vCon(e,:)) + I1x12;
    
    CMats.C1x21(mesh.tCon(e,:), mesh.vCon(e,:)) = ...
        CMats.C1x21(mesh.tCon(e,:), mesh.vCon(e,:)) + I1x21;
    
    CMats.C1x22(mesh.tCon(e,:), mesh.vCon(e,:)) = ...
        CMats.C1x22(mesh.tCon(e,:), mesh.vCon(e,:)) + I1x22;
    
    CMats.C2x11(mesh.tCon(e,:), mesh.vCon(e,:)) = ...
        CMats.C2x11(mesh.tCon(e,:), mesh.vCon(e,:)) + I2x11;
    
    CMats.C2x12(mesh.tCon(e,:), mesh.vCon(e,:)) = ...
        CMats.C2x12(mesh.tCon(e,:), mesh.vCon(e,:)) + I2x12;
    
    CMats.C2x21(mesh.tCon(e,:), mesh.vCon(e,:)) = ...
        CMats.C2x21(mesh.tCon(e,:), mesh.vCon(e,:)) + I2x21;
    
    CMats.C2x22(mesh.tCon(e,:), mesh.vCon(e,:)) = ...
        CMats.C2x22(mesh.tCon(e,:), mesh.vCon(e,:)) + I2x22;
    
end

end