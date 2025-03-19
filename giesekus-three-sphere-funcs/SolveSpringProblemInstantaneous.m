function [model,vels,press,stress] = SolveSpringProblemInstantaneous...
    (model,Xt,ft,FixedTime,FixedIts,isnewt,varargin)

%% setup

% unpack model parameters
ep = model.ep;

global FixedItsPrev calc;

if isnewt == 1
    calc = 1;
elseif isempty(FixedItsPrev) && isempty(varargin)
    calc = 1;
else
    calc = 0;
end

%% solve for finite element corrections

if isnewt == 1 || isempty(varargin)

    % if newtonian or t=0 do a Newtonian solve
    [U0,P0,T0] = SolveNewtonianProblem(model,FixedTime,FixedIts);

    % print
    fprintf('Solved Newtonian problem.\n');

else

    % we have information from previous timesteps and use it to intilise
    % the picard iterations
    U0 = varargin{1};
    P0 = varargin{2};
    T0 = varargin{3};
end

if isnewt == 0

    % problem is non-linear so run picard iterations
    if isempty(varargin)

        [U,P,T] = RunPicardIterations(model,U0,P0,T0,FixedTime,FixedIts);

    end

    % set solutions to be those at the final iteration
    U = U(:,end);
    P = P(:,end);
    T = T(:,end);

else

    % set final solutions to be the result of the Newtonian solve
    U = U0;
    P = P0;
    T = T0;
end

%% calculate full solution

% get mesh node coordinates
Xv   = [model.mesh.xv; model.mesh.yv];
Xp   = [model.mesh.xp; model.mesh.yp];
Xtau = [model.mesh.xt; model.mesh.yt];

% get uN at mesh nodes
for ii=1:3
    S        = RegStokesletVelocity(Xv,Xt(:,ii),ep);
    uN(:,ii) = S*ft(:,ii);
end
uN = sum(uN,2);

% get pN at mesh nodes
for ii=1:3
    PN  = RegStokesletPressure(Xp,Xt(:,ii),ep);
    pN(:,ii) = PN*ft(:,ii);
end
pN = sum(pN,2);

% get tN at mesh nodes
for ii=1:3
    dS    = RegStokesletDerivVelocity(Xtau,Xt(:,ii),ep);
    dS    = cell2mat(dS);
    dUNdx(:,ii) = dS(:,:,1)*ft(:,ii);
    dUNdy(:,ii) = dS(:,:,2)*ft(:,ii);

    [duNdx(:,ii),dvNdx(:,ii)] = ExtractVectorComponents(dUNdx(:,ii));
    [duNdy(:,ii),dvNdy(:,ii)] = ExtractVectorComponents(dUNdy(:,ii));

    tN11(:,ii) =  (duNdx(:,ii)+duNdx(:,ii));
    tN12(:,ii) = (duNdy(:,ii)+dvNdx(:,ii));
    tN21(:,ii) = (dvNdx(:,ii)+duNdy(:,ii));
    tN22(:,ii) =  (dvNdy(:,ii)+dvNdy(:,ii));

    tN(:,ii) = [tN11(:,ii); tN12(:,ii); tN21(:,ii); tN22(:,ii)];

end
tN = sum(tN,2);

% get final velocity and pressure at mesh nodes
u   = U + uN;
p   = P + pN;
tau = T + tN;

%% sort outputs

vels.u = u;
vels.U = U;
vels.uN = uN;

press.p = p;
press.P = P;
press.tN = pN;

stress.t = tau;
stress.T = T;
stress.tN = tN;

model.Xt = Xt;
model.ft = ft;
end