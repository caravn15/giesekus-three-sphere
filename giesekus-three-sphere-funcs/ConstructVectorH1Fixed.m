function H1FixedVecs = ConstructVectorH1Fixed(model,FixedIts)

% unpack required parameters
mesh  = model.mesh;
tN11  = FixedIts.tN11;
tN12  = FixedIts.tN12;
tN21  = FixedIts.tN21;
tN22  = FixedIts.tN22;
dUNdx = FixedIts.dUNdx;
dUNdy = FixedIts.dUNdy;
detJ  = FixedIts.detJ;
w2D   = FixedIts.w2D;
shape = FixedIts.shape;

% initilise full vectors
H1FixedVecs.H11 = sparse(model.mesh.nt,1);
H1FixedVecs.H12 = sparse(model.mesh.nt,1);
H1FixedVecs.H13 = sparse(model.mesh.nt,1);
H1FixedVecs.H14 = sparse(model.mesh.nt,1);

for ii=1:3

    % evaluate shape functions and stokeslet on each element
    for e=1:model.mesh.nelem

        % extract components
        [duNdx,dvNdx] = ExtractVectorComponents(dUNdx(:,e,ii));
        [duNdy,dvNdy] = ExtractVectorComponents(dUNdy(:,e,ii));

        % calculate integrals for each element
        [I1,I2,I3,I4] = ElementalStiffnessH1Fixed(shape(e,ii),...
            tN11(:,e,ii),tN12(:,e,ii),tN21(:,e,ii),tN22(:,e,ii),...
            duNdx,dvNdx,duNdy,dvNdy,detJ(:,e,ii),w2D);

        % assemble full vectors
        H1FixedVecs.H11(mesh.tCon(e,:), 1) = H1FixedVecs.H11...
            (mesh.tCon(e,:), 1) + I1;
        H1FixedVecs.H12(mesh.tCon(e,:), 1) = H1FixedVecs.H12...
            (mesh.tCon(e,:), 1) + I2;
        H1FixedVecs.H13(mesh.tCon(e,:), 1) = H1FixedVecs.H13...
            (mesh.tCon(e,:), 1) + I3;
        H1FixedVecs.H14(mesh.tCon(e,:), 1) = H1FixedVecs.H14...
            (mesh.tCon(e,:), 1) + I4;

    end
end

end