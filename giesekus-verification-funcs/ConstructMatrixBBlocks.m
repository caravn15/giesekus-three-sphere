function BMats = ConstructMatrixBBlocks(mesh,gdshape,shape,detJ,w2D)

% initilise full matrices
BMats.B1x11 = sparse(mesh.nv,mesh.nt);
BMats.B1x12 = sparse(mesh.nv,mesh.nt);
BMats.B1x21 = sparse(mesh.nv,mesh.nt);
BMats.B1x22 = sparse(mesh.nv,mesh.nt);

BMats.B2x11 = sparse(mesh.nv,mesh.nt);
BMats.B2x12 = sparse(mesh.nv,mesh.nt);
BMats.B2x21 = sparse(mesh.nv,mesh.nt);
BMats.B2x22 = sparse(mesh.nv,mesh.nt);

for e=1:mesh.nelem
    
    % calculate integrals for each element
    [I1x11,I1x12,I1x21,I1x22,I2x11,I2x12,I2x21,I2x22] = ...
        ElementalStiffnessB(gdshape{e},shape,detJ{e},w2D);
    
    % assemble full matrix
    BMats.B1x11(mesh.vCon(e,:), mesh.tCon(e,:)) = ...
        BMats.B1x11(mesh.vCon(e,:), mesh.tCon(e,:)) + I1x11;
    
    BMats.B1x12(mesh.vCon(e,:), mesh.tCon(e,:)) = ...
        BMats.B1x12(mesh.vCon(e,:), mesh.tCon(e,:)) + I1x12;
    
    BMats.B1x21(mesh.vCon(e,:), mesh.tCon(e,:)) = ...
        BMats.B1x21(mesh.vCon(e,:), mesh.tCon(e,:)) + I1x21;
    
    BMats.B1x22(mesh.vCon(e,:), mesh.tCon(e,:)) = ...
        BMats.B1x22(mesh.vCon(e,:), mesh.tCon(e,:)) + I1x22;
    
    BMats.B2x11(mesh.vCon(e,:), mesh.tCon(e,:)) = ...
        BMats.B2x11(mesh.vCon(e,:), mesh.tCon(e,:)) + I2x11;
    
    BMats.B2x12(mesh.vCon(e,:), mesh.tCon(e,:)) = ...
        BMats.B2x12(mesh.vCon(e,:), mesh.tCon(e,:)) + I2x12;
    
    BMats.B2x21(mesh.vCon(e,:), mesh.tCon(e,:)) = ...
        BMats.B2x21(mesh.vCon(e,:), mesh.tCon(e,:)) + I2x21;
    
    BMats.B2x22(mesh.vCon(e,:), mesh.tCon(e,:)) = ...
        BMats.B2x22(mesh.vCon(e,:), mesh.tCon(e,:)) + I2x22;
    
end

end