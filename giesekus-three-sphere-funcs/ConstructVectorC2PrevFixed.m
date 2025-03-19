function C2PrevFixedVecs = ConstructVectorC2PrevFixed(model,TptimeQuad,...
    FixedTime,FixedItsPrev)

% unpack required parameters
mesh  = model.mesh;
tN11  = FixedItsPrev.tN11;
tN12  = FixedItsPrev.tN12;
tN21  = FixedItsPrev.tN21;
tN22  = FixedItsPrev.tN22;
detJ  = FixedItsPrev.detJ;
w2D   = FixedItsPrev.w2D;
shape = FixedItsPrev.shape;

% initilise full vectors
C2PrevFixedVecs.C21 = sparse(model.mesh.nt,1);
C2PrevFixedVecs.C22 = sparse(model.mesh.nt,1);
C2PrevFixedVecs.C23 = sparse(model.mesh.nt,1);
C2PrevFixedVecs.C24 = sparse(model.mesh.nt,1);

for ii=1:3

    % evaluate shape functions and stokeslet on each element
    for e=1:model.mesh.nelem

        % calculate integrals for each element
        [I1,I2,I3,I4] = ElementalStiffnessC2PrevFixed(shape(e,ii),...
            tN11(:,e,ii),tN12(:,e,ii),tN21(:,e,ii),tN22(:,e,ii),...
            detJ(:,e,ii),w2D,TptimeQuad{e},FixedTime);

        % assemble full vectors
        C2PrevFixedVecs.C21(mesh.tCon(e,:), 1) = C2PrevFixedVecs.C21(...
            mesh.tCon(e,:), 1) + I1;
        C2PrevFixedVecs.C22(mesh.tCon(e,:), 1) = C2PrevFixedVecs.C22(...
            mesh.tCon(e,:), 1) + I2;
        C2PrevFixedVecs.C23(mesh.tCon(e,:), 1) = C2PrevFixedVecs.C23(...
            mesh.tCon(e,:), 1) + I3;
        C2PrevFixedVecs.C24(mesh.tCon(e,:), 1) = C2PrevFixedVecs.C24(...
            mesh.tCon(e,:), 1) + I4;

    end
end

end