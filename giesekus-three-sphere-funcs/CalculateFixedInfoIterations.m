function FixedIts = CalculateFixedInfoIterations(model,Xt,ft,...
    Tptime,isnewt,FixedTime)

%% setup

% unpack model parameters
ep = model.ep;

global FixedItsPrev calc;

if isnewt == 1
    calc = 1;
elseif isempty(FixedItsPrev) && model.tmin == 0
    calc = 1;
else
    calc = 0;
end

%% calculate information required for integral evaluation (sinh transform)

% setup quadrature rules
[x2D,w2D] = GetGaussianQuad2D(model.num.sinhquad);

for jj=1:3

    % get transformed points from reference element on each element
    [s,t] = TransformFromReference(Xt(:,jj),x2D,model);

    % get transformed points from global coordinates on each element
    [x,y] = TransformFromGlobal(Xt(:,jj),x2D,model);

    % get determinant of Jacobians of transformation
    detJ(:,:,jj) = CalculateSinhJacobian(Xt(:,jj),x2D,model);

    for e=1:model.mesh.nelem

        S = [s(:,e); t(:,e)];
        X = [x(:,e); y(:,e)];

        % calculate shape function at transfomred GQs
        shape(e,jj) = CalculateRefShape(S);

        % calculate shape function derivatives at transformed GQs
        dshape(e,jj) = CalculateRefShapeDeriv(S);

        % calculate stokeslet functions at transformed GQs (velocity)
        S  = RegStokesletVelocity(X,Xt(:,jj),ep);
        UN(:,e,jj) = S*ft(:,jj);

        % calculate stokeslet functions at transformed GQs (velocity
        % derivatives and stress)
        duN   = RegStokesletDerivVelocity(X,Xt(:,jj),ep);
        duN   = cell2mat(duN);
        dUNdx(:,e,jj) = duN(:,:,1)*ft(:,jj);
        dUNdy(:,e,jj) = duN(:,:,2)*ft(:,jj);

        [duNdx,dvNdx] = ExtractVectorComponents(dUNdx(:,e,jj));
        [duNdy,dvNdy] = ExtractVectorComponents(dUNdy(:,e,jj));

        tN11(:,e,jj) = duNdx+duNdx;
        tN12(:,e,jj) = duNdy+dvNdx;
        tN21(:,e,jj) = dvNdx+duNdy;
        tN22(:,e,jj) = dvNdy+dvNdy;

        % calculate stokeslet functions at transformed GQs (velocity
        % second derivatives and stress derivatives)
        d2uN = RegStokesletSecondDerivVelocity(X,Xt(:,jj),ep);
        d2uN   = cell2mat(d2uN);
        d2UNdxdx(:,e,jj) = d2uN(:,:,1,1)*ft(:,jj);
        d2UNdxdy(:,e,jj) = d2uN(:,:,1,2)*ft(:,jj);
        d2UNdydx(:,e,jj) = d2uN(:,:,2,1)*ft(:,jj);
        d2UNdydy(:,e,jj) = d2uN(:,:,2,2)*ft(:,jj);

        [d2uNdxdx,d2vNdxdx] = ExtractVectorComponents(d2UNdxdx(:,e,jj));
        [d2uNdxdy,d2vNdxdy] = ExtractVectorComponents(d2UNdxdy(:,e,jj));
        [d2uNdydx,d2vNdydx] = ExtractVectorComponents(d2UNdydx(:,e,jj));
        [d2uNdydy,d2vNdydy] = ExtractVectorComponents(d2UNdydy(:,e,jj));

        dtN11dx(:,e,jj) = d2uNdxdx + d2uNdxdx;
        dtN12dx(:,e,jj) = d2uNdxdy + d2vNdxdx;
        dtN21dx(:,e,jj) = d2vNdxdx + d2uNdxdy;
        dtN22dx(:,e,jj) = d2vNdxdy + d2vNdxdy;
        dtN11dy(:,e,jj) = d2uNdydx + d2uNdydx;
        dtN12dy(:,e,jj) = d2uNdydy + d2vNdydx;
        dtN21dy(:,e,jj) = d2vNdydx + d2uNdydy;
        dtN22dy(:,e,jj) = d2vNdydy + d2vNdydy;

    end
end

FixedIts.w2D     = w2D;
FixedIts.detJ    = detJ;
FixedIts.dshape  = dshape;
FixedIts.shape   = shape;
FixedIts.UN      = UN;
FixedIts.dUNdx   = dUNdx;
FixedIts.dUNdy   = dUNdy;
FixedIts.tN11    = tN11;
FixedIts.tN12    = tN12;
FixedIts.tN21    = tN21;
FixedIts.tN22    = tN22;
FixedIts.dtN11dx = dtN11dx;
FixedIts.dtN12dx = dtN12dx;
FixedIts.dtN21dx = dtN21dx;
FixedIts.dtN22dx = dtN22dx;
FixedIts.dtN11dy = dtN11dy;
FixedIts.dtN12dy = dtN12dy;
FixedIts.dtN21dy = dtN21dy;
FixedIts.dtN22dy = dtN22dy;

%% construct rhs vector components

% DVecs gives tau_N^m:D(v) - nabla(u_N^m):nabla(v)
DVecs = ConstructVectorD(model,FixedIts);

% package into structure
FixedIts.DVecs = DVecs;

% construct vector components that are only needed in the viscoelastic case
if isnewt == 0

    % C1FixedVecs gives (tau_N^m.tau_N^m):theta
    C1FixedVecs = ConstructVectorC1Fixed(model,FixedIts);

    % H1FixedVecs gives (tau_N^m.nabla(u_N^m)+(nabla(u_N^m))^T.tau_N^m):theta
    H1FixedVecs = ConstructVectorH1Fixed(model,FixedIts);

    % package into structure
    FixedIts.C1FixedVecs = C1FixedVecs;
    FixedIts.H1FixedVecs = H1FixedVecs;

    if calc == 0

        % evaluate stress from previous time at quadrature points
        TptimeQuad  = GetStressAtQuadPoints(model.mesh,...
            FixedTime.shape,Tptime);

        % C2CurFixedVecs gives tau_N^m:theta
        C2CurFixedVecs = ConstructVectorC2CurFixed(model,FixedIts);

        % C2PrevFixedVecs gives T^(m-1):theta+tau_N^(m-1):theta
        C2PrevFixedVecs = ConstructVectorC2PrevFixed(model,TptimeQuad,...
            FixedTime,FixedItsPrev);

        % H0FixedVecs gives (u_N^m.nabla(tau_N^m)):theta
        H0FixedVecs = ConstructVectorH0Fixed(model,FixedIts);

        % package into structure
        FixedIts.C2CurFixedVecs  = C2CurFixedVecs;
        FixedIts.C2PrevFixedVecs = C2PrevFixedVecs;
        FixedIts.H0FixedVecs     = H0FixedVecs;

    end

end

%% construct matrix components (only if viscoelastic)

if isnewt == 0

    % C1FixedMats gives (tau_N^m.T^(m,r)):theta
    C1FixedMats = ConstructMatrixC1FixedBlocks(model,FixedIts);

    % H1FixedMats gives (tau_N^m.nabla(U^(m,r))+(nabla(U^(m,r)))^T.tau_N^m):theta
    H1FixedMats = ConstructMatrixH1FixedBlocks(model,FixedIts);

    % package into structure
    FixedIts.C1FixedMats = C1FixedMats;
    FixedIts.H1FixedMats = H1FixedMats;

    if calc == 0

        % C2CurFixedMats gives T^(m,r):theta (already have this in C0Mats)
        C2CurFixedMats         = struct();
        C2CurFixedMats.C211x11 = FixedTime.C0Mats.C011x11;
        C2CurFixedMats.C211x12 = FixedTime.C0Mats.C011x12;
        C2CurFixedMats.C211x21 = FixedTime.C0Mats.C011x21;
        C2CurFixedMats.C211x22 = FixedTime.C0Mats.C011x22;
        C2CurFixedMats.C212x11 = FixedTime.C0Mats.C012x11;
        C2CurFixedMats.C212x12 = FixedTime.C0Mats.C012x12;
        C2CurFixedMats.C212x21 = FixedTime.C0Mats.C012x21;
        C2CurFixedMats.C212x22 = FixedTime.C0Mats.C012x22;
        C2CurFixedMats.C221x11 = FixedTime.C0Mats.C021x11;
        C2CurFixedMats.C221x12 = FixedTime.C0Mats.C021x12;
        C2CurFixedMats.C221x21 = FixedTime.C0Mats.C021x21;
        C2CurFixedMats.C221x22 = FixedTime.C0Mats.C021x22;
        C2CurFixedMats.C222x11 = FixedTime.C0Mats.C022x11;
        C2CurFixedMats.C222x12 = FixedTime.C0Mats.C022x12;
        C2CurFixedMats.C222x21 = FixedTime.C0Mats.C022x21;
        C2CurFixedMats.C222x22 = FixedTime.C0Mats.C022x22;

        % H0FixedMats gives U^(m,r).nabla(tau_N^m)
        H0FixedMats = ConstructMatrixH0FixedBlocks(model,FixedIts);

        % package into structure
        FixedIts.C2CurFixedMats = C2CurFixedMats;
        FixedIts.H0FixedMats    = H0FixedMats;

    end

end

%% boundary conditions

% calculate boundary conditions
bc = GetBoundaryConditionsStokeslet(model,Xt,ft);

% package into structure
FixedIts.bc = bc;

end