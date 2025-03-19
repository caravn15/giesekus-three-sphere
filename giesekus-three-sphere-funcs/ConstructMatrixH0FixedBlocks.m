function H0FixedMats = ConstructMatrixH0FixedBlocks(model,FixedIts)

% unpack required parameters
mesh    = model.mesh;
dtN11dx = FixedIts.dtN11dx;
dtN12dx = FixedIts.dtN12dx;
dtN21dx = FixedIts.dtN21dx;
dtN22dx = FixedIts.dtN22dx;
dtN11dy = FixedIts.dtN11dy;
dtN12dy = FixedIts.dtN12dy;
dtN21dy = FixedIts.dtN21dy;
dtN22dy = FixedIts.dtN22dy;
detJ    = FixedIts.detJ;
w2D     = FixedIts.w2D;
shape   = FixedIts.shape;

% initilise full matrices
H0FixedMats.H011x1 = sparse(mesh.nt,mesh.nv);
H0FixedMats.H012x1 = sparse(mesh.nt,mesh.nv);
H0FixedMats.H021x1 = sparse(mesh.nt,mesh.nv);
H0FixedMats.H022x1 = sparse(mesh.nt,mesh.nv);

H0FixedMats.H011x2 = sparse(mesh.nt,mesh.nv);
H0FixedMats.H012x2 = sparse(mesh.nt,mesh.nv);
H0FixedMats.H021x2 = sparse(mesh.nt,mesh.nv);
H0FixedMats.H022x2 = sparse(mesh.nt,mesh.nv);

for ii=1:3

    for e=1:mesh.nelem

        % calculate integrals for each element
        [I11x1,I12x1,I21x1,I22x1,I11x2,I12x2,I21x2,I22x2] = ...
            ElementalStiffnessH0FixedMat(shape(e,ii),...
            dtN11dx(:,e,ii),dtN12dx(:,e,ii),dtN21dx(:,e,ii),dtN22dx(:,e,ii),...
            dtN11dy(:,e,ii),dtN12dy(:,e,ii),dtN21dy(:,e,ii),dtN22dy(:,e,ii),...
            detJ(:,e,ii),w2D);

        % assemble full matrix
        H0FixedMats.H011x1(mesh.tCon(e,:), mesh.vCon(e,:)) = ...
            H0FixedMats.H011x1(mesh.tCon(e,:), mesh.vCon(e,:)) + I11x1;

        H0FixedMats.H012x1(mesh.tCon(e,:), mesh.vCon(e,:)) = ...
            H0FixedMats.H012x1(mesh.tCon(e,:), mesh.vCon(e,:)) + I12x1;

        H0FixedMats.H021x1(mesh.tCon(e,:), mesh.vCon(e,:)) = ...
            H0FixedMats.H021x1(mesh.tCon(e,:), mesh.vCon(e,:)) + I21x1;

        H0FixedMats.H022x1(mesh.tCon(e,:), mesh.vCon(e,:)) = ...
            H0FixedMats.H022x1(mesh.tCon(e,:), mesh.vCon(e,:)) + I22x1;

        H0FixedMats.H011x2(mesh.tCon(e,:), mesh.vCon(e,:)) = ...
            H0FixedMats.H011x2(mesh.tCon(e,:), mesh.vCon(e,:)) + I11x2;

        H0FixedMats.H012x2(mesh.tCon(e,:), mesh.vCon(e,:)) = ...
            H0FixedMats.H012x2(mesh.tCon(e,:), mesh.vCon(e,:)) + I12x2;

        H0FixedMats.H021x2(mesh.tCon(e,:), mesh.vCon(e,:)) = ...
            H0FixedMats.H021x2(mesh.tCon(e,:), mesh.vCon(e,:)) + I21x2;

        H0FixedMats.H022x2(mesh.tCon(e,:), mesh.vCon(e,:)) = ...
            H0FixedMats.H022x2(mesh.tCon(e,:), mesh.vCon(e,:)) + I22x2;

    end
end

end