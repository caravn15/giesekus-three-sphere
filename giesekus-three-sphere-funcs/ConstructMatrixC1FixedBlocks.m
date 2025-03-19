function C1FixedMats = ConstructMatrixC1FixedBlocks(model,FixedIts)

% unpack required parameters
mesh   = model.mesh;
tN11   = FixedIts.tN11;
tN12   = FixedIts.tN12;
tN21   = FixedIts.tN21;
tN22   = FixedIts.tN22;
detJ   = FixedIts.detJ;
w2D    = FixedIts.w2D;
shape  = FixedIts.shape;

% initilise full matrices
C1FixedMats.C111x11 = sparse(mesh.nt,mesh.nt);
C1FixedMats.C111x12 = sparse(mesh.nt,mesh.nt);
C1FixedMats.C111x21 = sparse(mesh.nt,mesh.nt);
C1FixedMats.C111x22 = sparse(mesh.nt,mesh.nt);

C1FixedMats.C112x11 = sparse(mesh.nt,mesh.nt);
C1FixedMats.C112x12 = sparse(mesh.nt,mesh.nt);
C1FixedMats.C112x21 = sparse(mesh.nt,mesh.nt);
C1FixedMats.C112x22 = sparse(mesh.nt,mesh.nt);

C1FixedMats.C121x11 = sparse(mesh.nt,mesh.nt);
C1FixedMats.C121x12 = sparse(mesh.nt,mesh.nt);
C1FixedMats.C121x21 = sparse(mesh.nt,mesh.nt);
C1FixedMats.C121x22 = sparse(mesh.nt,mesh.nt);

C1FixedMats.C122x11 = sparse(mesh.nt,mesh.nt);
C1FixedMats.C122x12 = sparse(mesh.nt,mesh.nt);
C1FixedMats.C122x21 = sparse(mesh.nt,mesh.nt);
C1FixedMats.C122x22 = sparse(mesh.nt,mesh.nt);

for ii=1:3

    for e=1:mesh.nelem

        % calculate integrals for each element
        [I11x11,I11x12,I11x21,I11x22,I12x11,I12x12,I12x21,I12x22,...
            I21x11,I21x12,I21x21,I21x22,I22x11,I22x12,I22x21,I22x22] = ...
            ElementalStiffnessC1FixedMat(shape(e,ii),tN11(:,e,ii),tN12(:,e,ii),...
            tN21(:,e,ii),tN22(:,e,ii),detJ(:,e,ii),w2D);

        % assemble full matrix
        C1FixedMats.C111x11(mesh.tCon(e,:), mesh.tCon(e,:)) = ...
            C1FixedMats.C111x11(mesh.tCon(e,:), mesh.tCon(e,:)) + I11x11;

        C1FixedMats.C111x12(mesh.tCon(e,:), mesh.tCon(e,:)) = ...
            C1FixedMats.C111x12(mesh.tCon(e,:), mesh.tCon(e,:)) + I11x12;

        C1FixedMats.C111x21(mesh.tCon(e,:), mesh.tCon(e,:)) = ...
            C1FixedMats.C111x21(mesh.tCon(e,:), mesh.tCon(e,:)) + I11x21;

        C1FixedMats.C111x22(mesh.tCon(e,:), mesh.tCon(e,:)) = ...
            C1FixedMats.C111x22(mesh.tCon(e,:), mesh.tCon(e,:)) + I11x22;

        C1FixedMats.C112x11(mesh.tCon(e,:), mesh.tCon(e,:)) = ...
            C1FixedMats.C112x11(mesh.tCon(e,:), mesh.tCon(e,:)) + I12x11;

        C1FixedMats.C112x12(mesh.tCon(e,:), mesh.tCon(e,:)) = ...
            C1FixedMats.C112x12(mesh.tCon(e,:), mesh.tCon(e,:)) + I12x12;

        C1FixedMats.C112x21(mesh.tCon(e,:), mesh.tCon(e,:)) = ...
            C1FixedMats.C112x21(mesh.tCon(e,:), mesh.tCon(e,:)) + I12x21;

        C1FixedMats.C112x22(mesh.tCon(e,:), mesh.tCon(e,:)) = ...
            C1FixedMats.C112x22(mesh.tCon(e,:), mesh.tCon(e,:)) + I12x22;

        C1FixedMats.C121x11(mesh.tCon(e,:), mesh.tCon(e,:)) = ...
            C1FixedMats.C121x11(mesh.tCon(e,:), mesh.tCon(e,:)) + I21x11;

        C1FixedMats.C121x12(mesh.tCon(e,:), mesh.tCon(e,:)) = ...
            C1FixedMats.C121x12(mesh.tCon(e,:), mesh.tCon(e,:)) + I21x12;

        C1FixedMats.C121x21(mesh.tCon(e,:), mesh.tCon(e,:)) = ...
            C1FixedMats.C121x21(mesh.tCon(e,:), mesh.tCon(e,:)) + I21x21;

        C1FixedMats.C121x22(mesh.tCon(e,:), mesh.tCon(e,:)) = ...
            C1FixedMats.C121x22(mesh.tCon(e,:), mesh.tCon(e,:)) + I21x22;

        C1FixedMats.C122x11(mesh.tCon(e,:), mesh.tCon(e,:)) = ...
            C1FixedMats.C122x11(mesh.tCon(e,:), mesh.tCon(e,:)) + I22x11;

        C1FixedMats.C122x12(mesh.tCon(e,:), mesh.tCon(e,:)) = ...
            C1FixedMats.C122x12(mesh.tCon(e,:), mesh.tCon(e,:)) + I22x12;

        C1FixedMats.C122x21(mesh.tCon(e,:), mesh.tCon(e,:)) = ...
            C1FixedMats.C122x21(mesh.tCon(e,:), mesh.tCon(e,:)) + I22x21;

        C1FixedMats.C122x22(mesh.tCon(e,:), mesh.tCon(e,:)) = ...
            C1FixedMats.C122x22(mesh.tCon(e,:), mesh.tCon(e,:)) + I22x22;

    end
end

end