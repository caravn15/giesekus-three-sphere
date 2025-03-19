function H0FixedVecs = ConstructVectorH0Fixed(model,FixedIts)

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
UN      = FixedIts.UN;
detJ    = FixedIts.detJ;
w2D     = FixedIts.w2D;
shape   = FixedIts.shape;

% initilise full vectors
H0FixedVecs.H01 = sparse(model.mesh.nt,1);
H0FixedVecs.H02 = sparse(model.mesh.nt,1);
H0FixedVecs.H03 = sparse(model.mesh.nt,1);
H0FixedVecs.H04 = sparse(model.mesh.nt,1);

for ii=1:3

    % evaluate shape functions and stokeslet on each element
    for e=1:model.mesh.nelem

        % extract components
        [uN,vN] = ExtractVectorComponents(UN(:,e,ii));

        % calculate integrals for each element
        [I1,I2,I3,I4] = ElementalStiffnessH0Fixed(shape(e,ii),...
            dtN11dx(:,e,ii),dtN12dx(:,e,ii),dtN21dx(:,e,ii),dtN22dx(:,e,ii),...
            dtN11dy(:,e,ii),dtN12dy(:,e,ii),dtN21dy(:,e,ii),dtN22dy(:,e,ii),...
            uN,vN,detJ(:,e,ii),w2D);

        % assemble full vectors
        H0FixedVecs.H01(mesh.tCon(e,:), 1) = H0FixedVecs.H01...
            (mesh.tCon(e,:), 1) + I1;
        H0FixedVecs.H02(mesh.tCon(e,:), 1) = H0FixedVecs.H02...
            (mesh.tCon(e,:), 1) + I2;
        H0FixedVecs.H03(mesh.tCon(e,:), 1) = H0FixedVecs.H03...
            (mesh.tCon(e,:), 1) + I3;
        H0FixedVecs.H04(mesh.tCon(e,:), 1) = H0FixedVecs.H04...
            (mesh.tCon(e,:), 1) + I4;

    end
end

end