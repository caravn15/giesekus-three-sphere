function H0Vecs = ConstructVectorH0(model,FixedIts,TpIt)

% unpack required parameters
xlim   = model.mesh.xlim;
ylim   = model.mesh.ylim;
mesh   = model.mesh;
UN     = FixedIts.UN;
detJ   = FixedIts.detJ;
w2D    = FixedIts.w2D;
dshape = FixedIts.dshape;
shape  = FixedIts.shape;

% initilise full vectors
H0Vecs.H01 = sparse(model.mesh.nt,1);
H0Vecs.H02 = sparse(model.mesh.nt,1);
H0Vecs.H03 = sparse(model.mesh.nt,1);
H0Vecs.H04 = sparse(model.mesh.nt,1);

for ii=1:3

    % evaluate shape functions and stokeslet on each element
    for e=1:model.mesh.nelem

        % extract components
        [uN,vN] = ExtractVectorComponents(UN(:,e,ii));

        % calculate global shape function derivatives at quad points
        gdshape(1,:,:) = 2/(xlim(e,2)-xlim(e,1))*dshape(e,ii).v(1,:,:);
        gdshape(2,:,:) = 2/(ylim(e,2)-ylim(e,1))*dshape(e,ii).v(2,:,:);

        % evaluate stress from previous iteration at quadrature points
        dTpItQuad  = GetFemStressDerivatives(TpIt,mesh,gdshape,e);

        % calculate integrals for each element
        [I1,I2,I3,I4] = ElementalStiffnessH0(shape(e,ii),dTpItQuad,...
            uN,vN,detJ(:,e,ii),w2D);

        % assemble full vectors
        H0Vecs.H01(mesh.tCon(e,:), 1) = H0Vecs.H01...
            (mesh.tCon(e,:), 1) + I1;
        H0Vecs.H02(mesh.tCon(e,:), 1) = H0Vecs.H02...
            (mesh.tCon(e,:), 1) + I2;
        H0Vecs.H03(mesh.tCon(e,:), 1) = H0Vecs.H03...
            (mesh.tCon(e,:), 1) + I3;
        H0Vecs.H04(mesh.tCon(e,:), 1) = H0Vecs.H04...
            (mesh.tCon(e,:), 1) + I4;

    end
end

end