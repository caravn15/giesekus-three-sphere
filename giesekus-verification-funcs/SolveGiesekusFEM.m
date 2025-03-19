function [u,pp,tau,mesh,uAn] = SolveGiesekusFEM(dom,params,nquad,tol,varargin)

%% setup
fprintf('RUNNING GIESEKUS VERIFICATION\n');
fprintf('==========================================\n')

% generate mesh
mesh = GenerateMeshRect(dom);

% calculate global information
[w2D,shape,detJ,gdshape] = CalculateGlobalInfoIterations(mesh,nquad);

if isempty(varargin)

    % solve Newtonian problem for initial guess
    [u0,p0,tau0,uAn] = SolveNewtonianProblem(mesh,params,w2D,shape,...
        gdshape,detJ);
else
    u0   = varargin{1};
    p0   = varargin{2};
    tau0 = varargin{3};
    uAn  = varargin{4};
end

%% solve

De = params.De;
uit(:,1) = u0;
pit(:,1) = p0;
tit(:,1) = tau0;

for ii=1:length(De)

    clear uit pit tit;

    params.De = De(ii);

    if ii==1
        uit(:,1) = u0;
        pit(:,1) = p0;
        tit(:,1) = tau0;
    else
        uit(:,1) = u;
        pit(:,1) = pp;
        tit(:,1) = tau;
    end

    r = 2;  epp = 1;
    while 1

        % calculate stress and velocity derivatives at quad points from
        % previous iteration
        TQuad  = GetStressAtQuadPoints(mesh,shape,tit(:,r-1));
        dUQuad = GetDerivVelocityAtQuadPoints(mesh,gdshape,uit(:,r-1));

        % calculate Jacobian and Residual vector
        J = CalculateNewtonsJacobian(mesh,params,dUQuad,TQuad,w2D,shape,...
            gdshape,detJ);
        [a,uAn] = CalculateNewtonsResidual(mesh,params);

        % solve linear system
        d = J\a;

        % extract velocity, pressure and stress
        uit(:,r) = d(1:2*mesh.nv);
        pit(:,r) = d(2*mesh.nv+1:2*mesh.nv+mesh.np);
        tit(:,r) = d(2*mesh.nv+mesh.np+1:end);

        if r>3
            epp = RunIterations(uit,tit,pit,r,epp);
        end

        uit(:,r) = epp*uit(:,r)+(1-epp)*uit(:,r-1);
        pit(:,r) = epp*pit(:,r)+(1-epp)*pit(:,r-1);
        tit(:,r) = epp*tit(:,r)+(1-epp)*tit(:,r-1);

        % calculate errors
        errv = norm(uit(:,r)-uit(:,r-1));
        errt = norm(tit(:,r)-tit(:,r-1));

        fprintf(['Iteration ',num2str(r-1),' errv = ',num2str(errv), ...
            ', errtau = ',num2str(errt),' De = ',num2str(params.De),'\n']);

        if errv < tol && errt < tol
            break;
        end
        r = r + 1;

    end

    u   = uit(:,end);
    pp   = pit(:,end);
    tau = tit(:,end);



end


end