function H1FixedMats = ConstructMatrixH1FixedBlocks(model,FixedIts)

% unpack required parameters
xlim   = model.mesh.xlim;
ylim   = model.mesh.ylim;
mesh   = model.mesh;
tN11   = FixedIts.tN11;
tN12   = FixedIts.tN12;
tN21   = FixedIts.tN21;
tN22   = FixedIts.tN22;
detJ   = FixedIts.detJ;
w2D    = FixedIts.w2D;
shape  = FixedIts.shape;
dshape = FixedIts.dshape;

% initilise full matrices
H1FixedMats.H111x1 = sparse(mesh.nt,mesh.nv);
H1FixedMats.H112x1 = sparse(mesh.nt,mesh.nv);
H1FixedMats.H121x1 = sparse(mesh.nt,mesh.nv);
H1FixedMats.H122x1 = sparse(mesh.nt,mesh.nv);

H1FixedMats.H111x2 = sparse(mesh.nt,mesh.nv);
H1FixedMats.H112x2 = sparse(mesh.nt,mesh.nv);
H1FixedMats.H121x2 = sparse(mesh.nt,mesh.nv);
H1FixedMats.H122x2 = sparse(mesh.nt,mesh.nv);

for ii=1:3

    for e=1:mesh.nelem

        % calculate global shape function derivatives at quad points
        gdshape(1,:,:) = 2/(xlim(e,2)-xlim(e,1))*dshape(e,ii).v(1,:,:);
        gdshape(2,:,:) = 2/(ylim(e,2)-ylim(e,1))*dshape(e,ii).v(2,:,:);

        % calculate integrals for each element
        [I11x1,I12x1,I21x1,I1x22,I11x2,I12x2,I21x2,I22x2] = ...
            ElementalStiffnessH1FixedMat(gdshape,shape(e,ii),tN11(:,e,ii),...
            tN12(:,e,ii),tN21(:,e,ii),tN22(:,e,ii),...
            detJ(:,e,ii),w2D);

        % assemble full matrix
        H1FixedMats.H111x1(mesh.tCon(e,:), mesh.vCon(e,:)) = ...
            H1FixedMats.H111x1(mesh.tCon(e,:), mesh.vCon(e,:)) + I11x1;

        H1FixedMats.H112x1(mesh.tCon(e,:), mesh.vCon(e,:)) = ...
            H1FixedMats.H112x1(mesh.tCon(e,:), mesh.vCon(e,:)) + I12x1;

        H1FixedMats.H121x1(mesh.tCon(e,:), mesh.vCon(e,:)) = ...
            H1FixedMats.H121x1(mesh.tCon(e,:), mesh.vCon(e,:)) + I21x1;

        H1FixedMats.H122x1(mesh.tCon(e,:), mesh.vCon(e,:)) = ...
            H1FixedMats.H122x1(mesh.tCon(e,:), mesh.vCon(e,:)) + I1x22;

        H1FixedMats.H111x2(mesh.tCon(e,:), mesh.vCon(e,:)) = ...
            H1FixedMats.H111x2(mesh.tCon(e,:), mesh.vCon(e,:)) + I11x2;

        H1FixedMats.H112x2(mesh.tCon(e,:), mesh.vCon(e,:)) = ...
            H1FixedMats.H112x2(mesh.tCon(e,:), mesh.vCon(e,:)) + I12x2;

        H1FixedMats.H121x2(mesh.tCon(e,:), mesh.vCon(e,:)) = ...
            H1FixedMats.H121x2(mesh.tCon(e,:), mesh.vCon(e,:)) + I21x2;

        H1FixedMats.H122x2(mesh.tCon(e,:), mesh.vCon(e,:)) = ...
            H1FixedMats.H122x2(mesh.tCon(e,:), mesh.vCon(e,:)) + I22x2;

    end
end

end