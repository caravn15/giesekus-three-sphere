function [Xt,velX] = TimestepEuler(model,vels,Xt,ft)

% extract components
xv = model.mesh.xv;
yv = model.mesh.yv;
dt = model.dt;
ep = model.ep;

% get fem velocity solution components
[U,V] = ExtractVectorComponents(vels.U);

% get fem velocity at point
InterpU = scatteredInterpolant(xv,yv,U);
InterpV = scatteredInterpolant(xv,yv,V);

for ii=1:3

    Ux0 = InterpU(Xt(1,ii),Xt(2,ii));
    Vx0 = InterpV(Xt(1,ii),Xt(2,ii));


    for jj=1:3
        % get stokeslet velocity at point
        S       = RegStokesletVelocity(Xt(:,ii),Xt(:,jj),ep);
        uN      = S*ft(:,jj);
        [u(:,jj),v(:,jj)]   = ExtractVectorComponents(uN);
    end

    u = sum(u,2);
    v = sum(v,2);

    % get total velocity at point
    u = Ux0 + u;
    v = Vx0 + v;
    U = [u;v];

    % time step point
    Xt(:,ii) = Xt(:,ii) + dt*U;
    velX(1,ii) = u;
    velX(2,ii) = v;
end

end