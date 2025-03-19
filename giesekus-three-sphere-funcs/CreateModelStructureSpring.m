function model = CreateModelStructureSpring(X0,kappa,Amp,chi,lr,ep,...
    alpha,Wi,dt,nt,tmax,tmin,nbeats,ngridvec)

model = struct( ...
    'dt'            , dt     , ... % timestep size
    'nt'            , nt     , ... % number of timesteps
    'tmax'          , tmax   , ... % maximum time to simulate to
    'tmin'          , tmin   , ... % start time
    'nbeats'        , nbeats , ... % number of beats
    'X0'            , X0     , ... % initial position of spheres
    'lr'            , lr     , ... % resting length of springs
    'ep'            , ep     , ... % regularisation parameter/sphere radii
    'kappa'         , kappa  , ... % ratio of spring stiffness to fluid stress
    'Amp'           , Amp    , ... % actuation amplitude
    'chi'           , chi    , ... % actuation phase difference
    'alpha'         , alpha  , ... % alpha parameter is Giesekus model
    'Wi'            , Wi     , ... % Weissenberg number in Giesekus model
    'mesh', ... 
    struct( ...
        'lx'        , [0,2]  , ... % domain dimensions in x
        'ly'        , [0,2]  , ... % domain dimensions in y
        'nx'        , ngridvec     , ... % # of elements in x direction
        'ny'        , ngridvec       ... % # of elements  in y direction
        ), ...
    'num', ...
    struct( ...
        'tol'       , 1e-3   , ... % tolerance for Picard iterations
        'gaussquad' , 5     , ... % # of quad points per direction (gauss) 15
        'sinhquad'  , 7       ... % # of quad points per direction (sinh) 15
        )...
    );
end
