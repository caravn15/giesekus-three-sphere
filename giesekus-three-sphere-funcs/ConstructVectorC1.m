function C1Vecs = ConstructVectorC1(model,FixedIts,TpIt)

% unpack required parameters
mesh  = model.mesh;
tN11  = FixedIts.tN11;
tN12  = FixedIts.tN12;
tN21  = FixedIts.tN21;
tN22  = FixedIts.tN22;
detJ  = FixedIts.detJ;
w2D   = FixedIts.w2D;
shape = FixedIts.shape;

% initilise full vectors
C1Vecs.C11 = sparse(model.mesh.nt,1);
C1Vecs.C12 = sparse(model.mesh.nt,1);
C1Vecs.C13 = sparse(model.mesh.nt,1);
C1Vecs.C14 = sparse(model.mesh.nt,1);

for ii=1:3

    % evaluate shape functions and stokeslet on each element
    for e=1:model.mesh.nelem

        % evaluate stress from previous iteration at quadrature points
        TpItQuad  = GetStressAtQuadPoints(mesh,shape(e,ii),TpIt,e);

        % calculate integrals for each element
        [I1,I2,I3,I4] = ElementalStiffnessC1(shape(e,ii),...
            tN11(:,e,ii),tN12(:,e,ii),tN21(:,e,ii),tN22(:,e,ii),...
            TpItQuad,detJ(:,e,ii),w2D);

        % assemble full vectors
        C1Vecs.C11(mesh.tCon(e,:), 1) = C1Vecs.C11(...
            mesh.tCon(e,:), 1) + I1;
        C1Vecs.C12(mesh.tCon(e,:), 1) = C1Vecs.C12(...
            mesh.tCon(e,:), 1) + I2;
        C1Vecs.C13(mesh.tCon(e,:), 1) = C1Vecs.C13(...
            mesh.tCon(e,:), 1) + I3;
        C1Vecs.C14(mesh.tCon(e,:), 1) = C1Vecs.C14(...
            mesh.tCon(e,:), 1) + I4;

    end

end