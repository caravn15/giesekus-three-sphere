function DVecs = ConstructVectorD(model,FixedIts)

% unpack required parameters
xlim   = model.mesh.xlim;
ylim   = model.mesh.ylim;
mesh   = model.mesh;
dshape = FixedIts.dshape;
dUNdx  = FixedIts.dUNdx;
dUNdy  = FixedIts.dUNdy;
tN11   = FixedIts.tN11;
tN12   = FixedIts.tN12;
tN21   = FixedIts.tN21;
tN22   = FixedIts.tN22;
detJ   = FixedIts.detJ;
w2D    = FixedIts.w2D;

% initilise full vectors
DVecs.D1 = sparse(model.mesh.nv,1);
DVecs.D2 = sparse(model.mesh.nv,1);

for ii=1:3

    % evaluate shape functions and stokeslet on each element
    for e=1:model.mesh.nelem

        % calculate global shape function derivatives at quad points
        gdshape(1,:,:) = 2/(xlim(e,2)-xlim(e,1))*dshape(e,ii).v(1,:,:);
        gdshape(2,:,:) = 2/(ylim(e,2)-ylim(e,1))*dshape(e,ii).v(2,:,:);

        % calculate integrals for each element
        [I1,I2] = ElementalStiffnessD(gdshape,dUNdx(:,e,ii),...
            dUNdy(:,e,ii),tN11(:,e,ii),tN12(:,e,ii),tN21(:,e,ii),...
            tN22(:,e,ii),detJ(:,e,ii),w2D);

        % assemble full vectors
        DVecs.D1(mesh.vCon(e,:), 1) = DVecs.D1(mesh.vCon(e,:), 1) + I1;
        DVecs.D2(mesh.vCon(e,:), 1) = DVecs.D2(mesh.vCon(e,:), 1) + I2;

    end
end

end