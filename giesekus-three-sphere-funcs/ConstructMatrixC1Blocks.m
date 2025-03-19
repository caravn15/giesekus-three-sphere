function C1Mats = ConstructMatrixC1Blocks(mesh,shape,TpItQuad,detJ,w2C1)

% initilise full matrices
C1Mats.C111x11 = sparse(mesh.nt,mesh.nt);
C1Mats.C111x12 = sparse(mesh.nt,mesh.nt);
C1Mats.C111x21 = sparse(mesh.nt,mesh.nt);
C1Mats.C111x22 = sparse(mesh.nt,mesh.nt);

C1Mats.C112x11 = sparse(mesh.nt,mesh.nt);
C1Mats.C112x12 = sparse(mesh.nt,mesh.nt);
C1Mats.C112x21 = sparse(mesh.nt,mesh.nt);
C1Mats.C112x22 = sparse(mesh.nt,mesh.nt);

C1Mats.C121x11 = sparse(mesh.nt,mesh.nt);
C1Mats.C121x12 = sparse(mesh.nt,mesh.nt);
C1Mats.C121x21 = sparse(mesh.nt,mesh.nt);
C1Mats.C121x22 = sparse(mesh.nt,mesh.nt);

C1Mats.C122x11 = sparse(mesh.nt,mesh.nt);
C1Mats.C122x12 = sparse(mesh.nt,mesh.nt);
C1Mats.C122x21 = sparse(mesh.nt,mesh.nt);
C1Mats.C122x22 = sparse(mesh.nt,mesh.nt);

for e=1:mesh.nelem
    
    % calculate integrals for each element
    [I11x11,I11x12,I11x21,I11x22,I12x11,I12x12,I12x21,I12x22,...
        I21x11,I21x12,I21x21,I21x22,I22x11,I22x12,I22x21,I22x22] = ...
        ElementalStiffnessC1Mat(TpItQuad{e},shape,detJ{e},w2C1);
    
    % assemble full matrix
    C1Mats.C111x11(mesh.tCon(e,:), mesh.tCon(e,:)) = ...
        C1Mats.C111x11(mesh.tCon(e,:), mesh.tCon(e,:)) + I11x11;
    
    C1Mats.C111x12(mesh.tCon(e,:), mesh.tCon(e,:)) = ...
        C1Mats.C111x12(mesh.tCon(e,:), mesh.tCon(e,:)) + I11x12;
    
    C1Mats.C111x21(mesh.tCon(e,:), mesh.tCon(e,:)) = ...
        C1Mats.C111x21(mesh.tCon(e,:), mesh.tCon(e,:)) + I11x21;
    
    C1Mats.C111x22(mesh.tCon(e,:), mesh.tCon(e,:)) = ...
        C1Mats.C111x22(mesh.tCon(e,:), mesh.tCon(e,:)) + I11x22;
    
    C1Mats.C112x11(mesh.tCon(e,:), mesh.tCon(e,:)) = ...
        C1Mats.C112x11(mesh.tCon(e,:), mesh.tCon(e,:)) + I12x11;
    
    C1Mats.C112x12(mesh.tCon(e,:), mesh.tCon(e,:)) = ...
        C1Mats.C112x12(mesh.tCon(e,:), mesh.tCon(e,:)) + I12x12;
    
    C1Mats.C112x21(mesh.tCon(e,:), mesh.tCon(e,:)) = ...
        C1Mats.C112x21(mesh.tCon(e,:), mesh.tCon(e,:)) + I12x21;
    
    C1Mats.C112x22(mesh.tCon(e,:), mesh.tCon(e,:)) = ...
        C1Mats.C112x22(mesh.tCon(e,:), mesh.tCon(e,:)) + I12x22;
    
    C1Mats.C121x11(mesh.tCon(e,:), mesh.tCon(e,:)) = ...
        C1Mats.C121x11(mesh.tCon(e,:), mesh.tCon(e,:)) + I21x11;
    
    C1Mats.C121x12(mesh.tCon(e,:), mesh.tCon(e,:)) = ...
        C1Mats.C121x12(mesh.tCon(e,:), mesh.tCon(e,:)) + I21x12;
    
    C1Mats.C121x21(mesh.tCon(e,:), mesh.tCon(e,:)) = ...
        C1Mats.C121x21(mesh.tCon(e,:), mesh.tCon(e,:)) + I21x21;
    
    C1Mats.C121x22(mesh.tCon(e,:), mesh.tCon(e,:)) = ...
        C1Mats.C121x22(mesh.tCon(e,:), mesh.tCon(e,:)) + I21x22;
    
    C1Mats.C122x11(mesh.tCon(e,:), mesh.tCon(e,:)) = ...
        C1Mats.C122x11(mesh.tCon(e,:), mesh.tCon(e,:)) + I22x11;
    
    C1Mats.C122x12(mesh.tCon(e,:), mesh.tCon(e,:)) = ...
        C1Mats.C122x12(mesh.tCon(e,:), mesh.tCon(e,:)) + I22x12;
    
    C1Mats.C122x21(mesh.tCon(e,:), mesh.tCon(e,:)) = ...
        C1Mats.C122x21(mesh.tCon(e,:), mesh.tCon(e,:)) + I22x21;
    
    C1Mats.C122x22(mesh.tCon(e,:), mesh.tCon(e,:)) = ...
        C1Mats.C122x22(mesh.tCon(e,:), mesh.tCon(e,:)) + I22x22;
    
end

end