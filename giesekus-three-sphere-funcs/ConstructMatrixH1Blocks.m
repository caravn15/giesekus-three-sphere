function H1Mats = ConstructMatrixH1Blocks(mesh,gdshape,shape,TpItQuad,detJ,w2D)

% initilise full matrices
H1Mats.H111x1 = sparse(mesh.nt,mesh.nv);
H1Mats.H112x1 = sparse(mesh.nt,mesh.nv);
H1Mats.H121x1 = sparse(mesh.nt,mesh.nv);
H1Mats.H122x1 = sparse(mesh.nt,mesh.nv);

H1Mats.H111x2 = sparse(mesh.nt,mesh.nv);
H1Mats.H112x2 = sparse(mesh.nt,mesh.nv);
H1Mats.H121x2 = sparse(mesh.nt,mesh.nv);
H1Mats.H122x2 = sparse(mesh.nt,mesh.nv);

for e=1:mesh.nelem
    
    % calculate integrals for each element
    [I11x1,I12x1,I21x1,I22x1,I11x2,I12x2,I21x2,I22x2] = ...
        ElementalStiffnessH1Mat(gdshape{e},shape,TpItQuad{e},detJ{e},w2D);
    
    % assemble full matrix
    H1Mats.H111x1(mesh.tCon(e,:), mesh.vCon(e,:)) = ...
        H1Mats.H111x1(mesh.tCon(e,:), mesh.vCon(e,:)) + I11x1;
    
    H1Mats.H112x1(mesh.tCon(e,:), mesh.vCon(e,:)) = ...
        H1Mats.H112x1(mesh.tCon(e,:), mesh.vCon(e,:)) + I12x1;
    
    H1Mats.H121x1(mesh.tCon(e,:), mesh.vCon(e,:)) = ...
        H1Mats.H121x1(mesh.tCon(e,:), mesh.vCon(e,:)) + I21x1;
    
    H1Mats.H122x1(mesh.tCon(e,:), mesh.vCon(e,:)) = ...
        H1Mats.H122x1(mesh.tCon(e,:), mesh.vCon(e,:)) + I22x1;
    
    H1Mats.H111x2(mesh.tCon(e,:), mesh.vCon(e,:)) = ...
        H1Mats.H111x2(mesh.tCon(e,:), mesh.vCon(e,:)) + I11x2;
    
    H1Mats.H112x2(mesh.tCon(e,:), mesh.vCon(e,:)) = ...
        H1Mats.H112x2(mesh.tCon(e,:), mesh.vCon(e,:)) + I12x2;
    
    H1Mats.H121x2(mesh.tCon(e,:), mesh.vCon(e,:)) = ...
        H1Mats.H121x2(mesh.tCon(e,:), mesh.vCon(e,:)) + I21x2;
    
    H1Mats.H122x2(mesh.tCon(e,:), mesh.vCon(e,:)) = ...
        H1Mats.H122x2(mesh.tCon(e,:), mesh.vCon(e,:)) + I22x2;
    
end

end