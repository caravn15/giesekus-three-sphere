function H1Vecs = ConstructVectorH1(model,FixedIts,TpIt)

% unpack required parameters
mesh  = model.mesh;
dUNdx = FixedIts.dUNdx;
dUNdy = FixedIts.dUNdy;
detJ  = FixedIts.detJ;
w2D   = FixedIts.w2D;
shape = FixedIts.shape;

% initilise full vectors
H1Vecs.H11 = sparse(model.mesh.nt,1);
H1Vecs.H12 = sparse(model.mesh.nt,1);
H1Vecs.H13 = sparse(model.mesh.nt,1);
H1Vecs.H14 = sparse(model.mesh.nt,1);

for ii=1:3

    % evaluate shape functions and stokeslet on each element
    for e=1:model.mesh.nelem

        % evaluate stress from previous iteration at quadrature points
        TpItQuad  = GetStressAtQuadPoints(mesh,shape(e,ii),TpIt,e);

        % extract components
        [duNdx,dvNdx] = ExtractVectorComponents(dUNdx(:,e,ii));
        [duNdy,dvNdy] = ExtractVectorComponents(dUNdy(:,e,ii));

        % calculate integrals for each element
        [I1,I2,I3,I4] = ElementalStiffnessH1(shape(e,ii),...
           TpItQuad,duNdx,dvNdx,duNdy,dvNdy,detJ(:,e,ii),w2D);

        % assemble full vectors
        H1Vecs.H11(mesh.tCon(e,:), 1) = H1Vecs.H11...
            (mesh.tCon(e,:), 1) + I1;
        H1Vecs.H12(mesh.tCon(e,:), 1) = H1Vecs.H12...
            (mesh.tCon(e,:), 1) + I2;
        H1Vecs.H13(mesh.tCon(e,:), 1) = H1Vecs.H13...
            (mesh.tCon(e,:), 1) + I3;
        H1Vecs.H14(mesh.tCon(e,:), 1) = H1Vecs.H14...
            (mesh.tCon(e,:), 1) + I4;

    end
end

end