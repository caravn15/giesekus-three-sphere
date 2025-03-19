function DMats = ConstructMatrixDBlocks(mesh,shape,TQuad,detJ,w2D)

% initilise full matrices
DMats.D11x11 = sparse(mesh.nt,mesh.nt);
DMats.D11x12 = sparse(mesh.nt,mesh.nt);
DMats.D11x21 = sparse(mesh.nt,mesh.nt);
DMats.D11x22 = sparse(mesh.nt,mesh.nt);

DMats.D12x11 = sparse(mesh.nt,mesh.nt);
DMats.D12x12 = sparse(mesh.nt,mesh.nt);
DMats.D12x21 = sparse(mesh.nt,mesh.nt);
DMats.D12x22 = sparse(mesh.nt,mesh.nt);

DMats.D21x11 = sparse(mesh.nt,mesh.nt);
DMats.D21x12 = sparse(mesh.nt,mesh.nt);
DMats.D21x21 = sparse(mesh.nt,mesh.nt);
DMats.D21x22 = sparse(mesh.nt,mesh.nt);

DMats.D22x11 = sparse(mesh.nt,mesh.nt);
DMats.D22x12 = sparse(mesh.nt,mesh.nt);
DMats.D22x21 = sparse(mesh.nt,mesh.nt);
DMats.D22x22 = sparse(mesh.nt,mesh.nt);

for e=1:mesh.nelem
    
    % calculate integrals for each element
    [I11x11,I11x12,I11x21,I11x22,I12x11,I12x12,I12x21,I12x22,...
        I21x11,I21x12,I21x21,I21x22,I22x11,I22x12,I22x21,I22x22] = ...
        ElementalStiffnessD(TQuad{e},shape,detJ{e},w2D);
    
    % assemble full matrix
    DMats.D11x11(mesh.tCon(e,:), mesh.tCon(e,:)) = ...
        DMats.D11x11(mesh.tCon(e,:), mesh.tCon(e,:)) + I11x11;
    
    DMats.D11x12(mesh.tCon(e,:), mesh.tCon(e,:)) = ...
        DMats.D11x12(mesh.tCon(e,:), mesh.tCon(e,:)) + I11x12;
    
    DMats.D11x21(mesh.tCon(e,:), mesh.tCon(e,:)) = ...
        DMats.D11x21(mesh.tCon(e,:), mesh.tCon(e,:)) + I11x21;
    
    DMats.D11x22(mesh.tCon(e,:), mesh.tCon(e,:)) = ...
        DMats.D11x22(mesh.tCon(e,:), mesh.tCon(e,:)) + I11x22;
    
    DMats.D12x11(mesh.tCon(e,:), mesh.tCon(e,:)) = ...
        DMats.D12x11(mesh.tCon(e,:), mesh.tCon(e,:)) + I12x11;
    
    DMats.D12x12(mesh.tCon(e,:), mesh.tCon(e,:)) = ...
        DMats.D12x12(mesh.tCon(e,:), mesh.tCon(e,:)) + I12x12;
    
    DMats.D12x21(mesh.tCon(e,:), mesh.tCon(e,:)) = ...
        DMats.D12x21(mesh.tCon(e,:), mesh.tCon(e,:)) + I12x21;
    
    DMats.D12x22(mesh.tCon(e,:), mesh.tCon(e,:)) = ...
        DMats.D12x22(mesh.tCon(e,:), mesh.tCon(e,:)) + I12x22;
    
    DMats.D21x11(mesh.tCon(e,:), mesh.tCon(e,:)) = ...
        DMats.D21x11(mesh.tCon(e,:), mesh.tCon(e,:)) + I21x11;
    
    DMats.D21x12(mesh.tCon(e,:), mesh.tCon(e,:)) = ...
        DMats.D21x12(mesh.tCon(e,:), mesh.tCon(e,:)) + I21x12;
    
    DMats.D21x21(mesh.tCon(e,:), mesh.tCon(e,:)) = ...
        DMats.D21x21(mesh.tCon(e,:), mesh.tCon(e,:)) + I21x21;
    
    DMats.D21x22(mesh.tCon(e,:), mesh.tCon(e,:)) = ...
        DMats.D21x22(mesh.tCon(e,:), mesh.tCon(e,:)) + I21x22;
    
    DMats.D22x11(mesh.tCon(e,:), mesh.tCon(e,:)) = ...
        DMats.D22x11(mesh.tCon(e,:), mesh.tCon(e,:)) + I22x11;
    
    DMats.D22x12(mesh.tCon(e,:), mesh.tCon(e,:)) = ...
        DMats.D22x12(mesh.tCon(e,:), mesh.tCon(e,:)) + I22x12;
    
    DMats.D22x21(mesh.tCon(e,:), mesh.tCon(e,:)) = ...
        DMats.D22x21(mesh.tCon(e,:), mesh.tCon(e,:)) + I22x21;
    
    DMats.D22x22(mesh.tCon(e,:), mesh.tCon(e,:)) = ...
        DMats.D22x22(mesh.tCon(e,:), mesh.tCon(e,:)) + I22x22;
    
end

end