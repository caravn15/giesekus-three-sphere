function C2CurFixedVecs = ConstructVectorC2CurFixed(model,FixedIts)

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
C2CurFixedVecs.C21 = sparse(model.mesh.nt,1);
C2CurFixedVecs.C22 = sparse(model.mesh.nt,1);
C2CurFixedVecs.C23 = sparse(model.mesh.nt,1);
C2CurFixedVecs.C24 = sparse(model.mesh.nt,1);

for ii=1:3

    % evaluate shape functions and stokeslet on each element
    for e=1:model.mesh.nelem

        % calculate integrals for each element
        [I1,I2,I3,I4] = ElementalStiffnessC2CurFixed(shape(e,ii),...
            tN11(:,e,ii),tN12(:,e,ii),tN21(:,e,ii),tN22(:,e,ii),...
            detJ(:,e,ii),w2D);

        % assemble full vectors
        C2CurFixedVecs.C21(mesh.tCon(e,:), 1) = C2CurFixedVecs.C21(...
            mesh.tCon(e,:), 1) + I1;
        C2CurFixedVecs.C22(mesh.tCon(e,:), 1) = C2CurFixedVecs.C22(...
            mesh.tCon(e,:), 1) + I2;
        C2CurFixedVecs.C23(mesh.tCon(e,:), 1) = C2CurFixedVecs.C23(...
            mesh.tCon(e,:), 1) + I3;
        C2CurFixedVecs.C24(mesh.tCon(e,:), 1) = C2CurFixedVecs.C24(...
            mesh.tCon(e,:), 1) + I4;

    end
end

end