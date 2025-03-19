function J = AssembleJacobianBlocks(mesh,params,AMats,BMats,CMats,DMats,FMats)

% extract parameters
De    = params.De;
alpha = params.alpha;

% zero blocks
Zvv = zeros(2*mesh.nv,2*mesh.nv);
Zpp = zeros(mesh.np,mesh.np);
Zpt = zeros(mesh.np,4*mesh.nt);

% A blocks
A  = [AMats.A1 AMats.A2];

% B blocks
B = [BMats.B1x11 BMats.B1x12 BMats.B1x21 BMats.B1x22 ;
     BMats.B2x11 BMats.B2x12 BMats.B2x21 BMats.B2x22];

% C blocks
C = [CMats.C1x11 CMats.C2x11 ;
     CMats.C1x12 CMats.C2x12 ;
     CMats.C1x21 CMats.C2x21 ;
     CMats.C1x22 CMats.C2x22];

% D blocks
D = [DMats.D11x11 DMats.D11x12 DMats.D11x21 DMats.D11x22 ;
     DMats.D12x11 DMats.D12x12 DMats.D12x21 DMats.D12x22 ;
     DMats.D21x11 DMats.D21x12 DMats.D21x21 DMats.D21x22 ;
     DMats.D22x11 DMats.D22x12 DMats.D22x21 DMats.D22x22];

% F blocks
F = [FMats.F11x11 FMats.F11x12 FMats.F11x21 FMats.F11x22 ;
     FMats.F12x11 FMats.F12x12 FMats.F12x21 FMats.F12x22 ;
     FMats.F21x11 FMats.F21x12 FMats.F21x21 FMats.F21x22 ;
     FMats.F22x11 FMats.F22x12 FMats.F22x21 FMats.F22x22];

% velocity block in Giesekus model
Gu = -2*B' - De*C;

% stress block in Giesekus model
Gt = alpha*De*D+F;

% assemble full matrix
J = [Zvv -A'   B    ;
     A    Zpp  Zpt  ;
     Gu   Zpt' Gt  ];

end