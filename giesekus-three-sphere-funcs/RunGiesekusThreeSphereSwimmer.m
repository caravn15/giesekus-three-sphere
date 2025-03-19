function [model,Xt,ft,velX,Utime,Ptime,Ttime,Uptime,Pptime,Tptime,FinalFixedIts,vels] = ...
    RunGiesekusThreeSphereSwimmer(X0,kappa,Amp,chi,lr,ep,alpha,Wi,dt,nbeats,tmin,ngridvec,varargin)

%% setup

global FixedItsPrev;

if tmin == 0
    FixedItsPrev = [];
end

tvec = tmin:dt:tmin+nbeats;
nt   = length(tvec);
tmax = tvec(end);

if nbeats == 0
    tmax = 0;
    nt = 1;
    tvec = 0;

end

% create model structure
model = CreateModelStructureSpring(X0,kappa,Amp,chi,lr,ep,alpha,Wi,...
    dt,nt,tmax,tmin,nbeats,ngridvec);

% create newtonian indicator
if model.Wi == 0 || model.alpha == 0
    isnewt = 1;
else
    isnewt = 0;
end

% generate mesh
model.mesh = GenerateMeshRect(model);

% print parameters
PrintParameters(model);

% calculate information that is unchanged over time
fprintf('Calculating blocks fixed across time...')
FixedTime = CalculateFixedInfoTime(model);
fprintf('Solved!\n');

%% solve for swimmer evolution

Xt(:,:,1) = model.X0;

for ii=1:nt

    fprintf(['Solving problem for timestep ',num2str(ii),'/',...,
        num2str(nt),':\n']);

    % calculate forces
    ft(:,:,ii) = CalculateForces(model,Xt(:,:,ii),tvec(ii));

    if ii == 1 && model.tmin == 0

        Tptime   = [];

        % calculate blocks fixed across iterations
        FixedIts = CalculateFixedInfoIterations(model,Xt(:,:,ii),...
            ft(:,:,ii),Tptime,isnewt,FixedTime);

        % solve problem
        [model,vels{ii},press{ii},stress{ii}] = ...
            SolveSpringProblemInstantaneous(model,Xt(:,:,ii),...
            ft(:,:,ii),FixedTime,FixedIts,...
            isnewt);

        Tptime = stress{ii}.T;
        Uptime  = vels{ii}.U;
        Pptime  = press{ii}.P;

    else

        if ~isempty(varargin) && ii == 1
            Tptime = varargin{1}.Tptime;
            Uptime = varargin{1}.Uptime;
            Pptime = varargin{1}.Pptime;

        end

        if ii ~=1
            FinalFixedIts = FixedIts;
        end

        FixedIts = CalculateFixedInfoIterations(model,Xt(:,:,ii),...
            ft(:,:,ii),Tptime,isnewt,FixedTime);

        % solve  problem
        [model,vels{ii},press{ii},stress{ii}] = ...
            SolveSpringProblemInstantaneous(model,Xt(:,:,ii),...
            ft(:,:,ii),FixedTime,FixedIts,...
            isnewt);


        Tptime = stress{ii}.T;
        Uptime  = vels{ii}.U;
        Pptime  = press{ii}.P;

    end

    % timestep sphere positions
    [Xt(:,:,ii+1),velX(:,:,ii)] = TimestepEuler(model,vels{ii},Xt(:,:,ii),...
        ft(:,:,ii));

    FixedItsPrev = FixedIts;

    fprintf('Solved!\n-----------------------------------------------\n');

end

Utime = [];
Ptime = [];
Ttime = [];

for ii=1:25:nt
    Utime = [Utime, vels{ii}.U];
    Ptime = [Ptime, press{ii}.P];
    Ttime = [Ttime, stress{ii}.T];
end

if ii~=1

    FinalFixedIts2.tN11 = FinalFixedIts.tN11;
    FinalFixedIts2.tN12 = FinalFixedIts.tN12;
    FinalFixedIts2.tN21 = FinalFixedIts.tN21;
    FinalFixedIts2.tN22 = FinalFixedIts.tN22;
    FinalFixedIts2.detJ = FinalFixedIts.detJ;
    FinalFixedIts2.w2D = FinalFixedIts.w2D;
    FinalFixedIts2.shape = FinalFixedIts.shape;
    model2.tmax = model.tmax;
    model2.mesh = model.mesh;

    clear FinalFixedIts vels1 velsend model;

    vels1 = [];
    velsend = [];
    FinalFixedIts = FinalFixedIts2;
    model = model2;

    clear FinalFixedIts2 model2;

    vels = [];

else
    FinalFixedIts = [];
end
end