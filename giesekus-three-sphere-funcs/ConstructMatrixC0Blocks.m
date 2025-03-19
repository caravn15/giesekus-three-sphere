function C0Mats = ConstructMatrixC0Blocks(mesh,shape,detJ,w2D)

% initilise full matrices
C0Mats.C011x11 = sparse(mesh.nt,mesh.nt);
C0Mats.C011x12 = sparse(mesh.nt,mesh.nt);
C0Mats.C011x21 = sparse(mesh.nt,mesh.nt);
C0Mats.C011x22 = sparse(mesh.nt,mesh.nt);

C0Mats.C012x11 = sparse(mesh.nt,mesh.nt);
C0Mats.C012x12 = sparse(mesh.nt,mesh.nt);
C0Mats.C012x21 = sparse(mesh.nt,mesh.nt);
C0Mats.C012x22 = sparse(mesh.nt,mesh.nt);

C0Mats.C021x11 = sparse(mesh.nt,mesh.nt);
C0Mats.C021x12 = sparse(mesh.nt,mesh.nt);
C0Mats.C021x21 = sparse(mesh.nt,mesh.nt);
C0Mats.C021x22 = sparse(mesh.nt,mesh.nt);

C0Mats.C022x11 = sparse(mesh.nt,mesh.nt);
C0Mats.C022x12 = sparse(mesh.nt,mesh.nt);
C0Mats.C022x21 = sparse(mesh.nt,mesh.nt);
C0Mats.C022x22 = sparse(mesh.nt,mesh.nt);

for e=1:mesh.nelem
    
    % calculate integrals for each element
    [I11x11,I11x12,I11x21,I11x22,I12x11,I12x12,I12x21,I12x22,...
        I21x11,I21x12,I21x21,I21x22,I22x11,I22x12,I22x21,I22x22] = ...
        ElementalStiffnessC0(shape,detJ{e},w2D);
    
    % assemble full matrix
    C0Mats.C011x11(mesh.tCon(e,:), mesh.tCon(e,:)) = ...
        C0Mats.C011x11(mesh.tCon(e,:), mesh.tCon(e,:)) + I11x11;
    
    C0Mats.C011x12(mesh.tCon(e,:), mesh.tCon(e,:)) = ...
        C0Mats.C011x12(mesh.tCon(e,:), mesh.tCon(e,:)) + I11x12;
    
    C0Mats.C011x21(mesh.tCon(e,:), mesh.tCon(e,:)) = ...
        C0Mats.C011x21(mesh.tCon(e,:), mesh.tCon(e,:)) + I11x21;
    
    C0Mats.C011x22(mesh.tCon(e,:), mesh.tCon(e,:)) = ...
        C0Mats.C011x22(mesh.tCon(e,:), mesh.tCon(e,:)) + I11x22;
    
    C0Mats.C012x11(mesh.tCon(e,:), mesh.tCon(e,:)) = ...
        C0Mats.C012x11(mesh.tCon(e,:), mesh.tCon(e,:)) + I12x11;
    
    C0Mats.C012x12(mesh.tCon(e,:), mesh.tCon(e,:)) = ...
        C0Mats.C012x12(mesh.tCon(e,:), mesh.tCon(e,:)) + I12x12;
    
    C0Mats.C012x21(mesh.tCon(e,:), mesh.tCon(e,:)) = ...
        C0Mats.C012x21(mesh.tCon(e,:), mesh.tCon(e,:)) + I12x21;
    
    C0Mats.C012x22(mesh.tCon(e,:), mesh.tCon(e,:)) = ...
        C0Mats.C012x22(mesh.tCon(e,:), mesh.tCon(e,:)) + I12x22;
    
    C0Mats.C021x11(mesh.tCon(e,:), mesh.tCon(e,:)) = ...
        C0Mats.C021x11(mesh.tCon(e,:), mesh.tCon(e,:)) + I21x11;
    
    C0Mats.C021x12(mesh.tCon(e,:), mesh.tCon(e,:)) = ...
        C0Mats.C021x12(mesh.tCon(e,:), mesh.tCon(e,:)) + I21x12;
    
    C0Mats.C021x21(mesh.tCon(e,:), mesh.tCon(e,:)) = ...
        C0Mats.C021x21(mesh.tCon(e,:), mesh.tCon(e,:)) + I21x21;
    
    C0Mats.C021x22(mesh.tCon(e,:), mesh.tCon(e,:)) = ...
        C0Mats.C021x22(mesh.tCon(e,:), mesh.tCon(e,:)) + I21x22;
    
    C0Mats.C022x11(mesh.tCon(e,:), mesh.tCon(e,:)) = ...
        C0Mats.C022x11(mesh.tCon(e,:), mesh.tCon(e,:)) + I22x11;
    
    C0Mats.C022x12(mesh.tCon(e,:), mesh.tCon(e,:)) = ...
        C0Mats.C022x12(mesh.tCon(e,:), mesh.tCon(e,:)) + I22x12;
    
    C0Mats.C022x21(mesh.tCon(e,:), mesh.tCon(e,:)) = ...
        C0Mats.C022x21(mesh.tCon(e,:), mesh.tCon(e,:)) + I22x21;
    
    C0Mats.C022x22(mesh.tCon(e,:), mesh.tCon(e,:)) = ...
        C0Mats.C022x22(mesh.tCon(e,:), mesh.tCon(e,:)) + I22x22;
    
end

end