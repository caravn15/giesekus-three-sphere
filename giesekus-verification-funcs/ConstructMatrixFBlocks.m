function FMats = ConstructMatrixFBlocks(mesh,shape,detJ,w2D)

% initilise full matrices
FMats.F11x11 = sparse(mesh.nt,mesh.nt);
FMats.F11x12 = sparse(mesh.nt,mesh.nt);
FMats.F11x21 = sparse(mesh.nt,mesh.nt);
FMats.F11x22 = sparse(mesh.nt,mesh.nt);

FMats.F12x11 = sparse(mesh.nt,mesh.nt);
FMats.F12x12 = sparse(mesh.nt,mesh.nt);
FMats.F12x21 = sparse(mesh.nt,mesh.nt);
FMats.F12x22 = sparse(mesh.nt,mesh.nt);

FMats.F21x11 = sparse(mesh.nt,mesh.nt);
FMats.F21x12 = sparse(mesh.nt,mesh.nt);
FMats.F21x21 = sparse(mesh.nt,mesh.nt);
FMats.F21x22 = sparse(mesh.nt,mesh.nt);

FMats.F22x11 = sparse(mesh.nt,mesh.nt);
FMats.F22x12 = sparse(mesh.nt,mesh.nt);
FMats.F22x21 = sparse(mesh.nt,mesh.nt);
FMats.F22x22 = sparse(mesh.nt,mesh.nt);

for e=1:mesh.nelem
    
    % calculate integrals for each element
    [I11x11,I11x12,I11x21,I11x22,I12x11,I12x12,I12x21,I12x22,...
        I21x11,I21x12,I21x21,I21x22,I22x11,I22x12,I22x21,I22x22] = ...
        ElementalStiffnessF(shape,detJ{e},w2D);
    
    % assemble full matrix
    FMats.F11x11(mesh.tCon(e,:), mesh.tCon(e,:)) = ...
        FMats.F11x11(mesh.tCon(e,:), mesh.tCon(e,:)) + I11x11;
    
    FMats.F11x12(mesh.tCon(e,:), mesh.tCon(e,:)) = ...
        FMats.F11x12(mesh.tCon(e,:), mesh.tCon(e,:)) + I11x12;
    
    FMats.F11x21(mesh.tCon(e,:), mesh.tCon(e,:)) = ...
        FMats.F11x21(mesh.tCon(e,:), mesh.tCon(e,:)) + I11x21;
    
    FMats.F11x22(mesh.tCon(e,:), mesh.tCon(e,:)) = ...
        FMats.F11x22(mesh.tCon(e,:), mesh.tCon(e,:)) + I11x22;
    
    FMats.F12x11(mesh.tCon(e,:), mesh.tCon(e,:)) = ...
        FMats.F12x11(mesh.tCon(e,:), mesh.tCon(e,:)) + I12x11;
    
    FMats.F12x12(mesh.tCon(e,:), mesh.tCon(e,:)) = ...
        FMats.F12x12(mesh.tCon(e,:), mesh.tCon(e,:)) + I12x12;
    
    FMats.F12x21(mesh.tCon(e,:), mesh.tCon(e,:)) = ...
        FMats.F12x21(mesh.tCon(e,:), mesh.tCon(e,:)) + I12x21;
    
    FMats.F12x22(mesh.tCon(e,:), mesh.tCon(e,:)) = ...
        FMats.F12x22(mesh.tCon(e,:), mesh.tCon(e,:)) + I12x22;
    
    FMats.F21x11(mesh.tCon(e,:), mesh.tCon(e,:)) = ...
        FMats.F21x11(mesh.tCon(e,:), mesh.tCon(e,:)) + I21x11;
    
    FMats.F21x12(mesh.tCon(e,:), mesh.tCon(e,:)) = ...
        FMats.F21x12(mesh.tCon(e,:), mesh.tCon(e,:)) + I21x12;
    
    FMats.F21x21(mesh.tCon(e,:), mesh.tCon(e,:)) = ...
        FMats.F21x21(mesh.tCon(e,:), mesh.tCon(e,:)) + I21x21;
    
    FMats.F21x22(mesh.tCon(e,:), mesh.tCon(e,:)) = ...
        FMats.F21x22(mesh.tCon(e,:), mesh.tCon(e,:)) + I21x22;
    
    FMats.F22x11(mesh.tCon(e,:), mesh.tCon(e,:)) = ...
        FMats.F22x11(mesh.tCon(e,:), mesh.tCon(e,:)) + I22x11;
    
    FMats.F22x12(mesh.tCon(e,:), mesh.tCon(e,:)) = ...
        FMats.F22x12(mesh.tCon(e,:), mesh.tCon(e,:)) + I22x12;
    
    FMats.F22x21(mesh.tCon(e,:), mesh.tCon(e,:)) = ...
        FMats.F22x21(mesh.tCon(e,:), mesh.tCon(e,:)) + I22x21;
    
    FMats.F22x22(mesh.tCon(e,:), mesh.tCon(e,:)) = ...
        FMats.F22x22(mesh.tCon(e,:), mesh.tCon(e,:)) + I22x22;
    
end

end