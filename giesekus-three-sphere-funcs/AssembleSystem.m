function [J,a] = AssembleSystem(alpha,Wi,dt,model,isunsteady,H0M,H1M,...
    BM,AM,C0M,C1M,C2M,C1V,C2V,H0V,H1V,DV)

% zero blocks
Zvv = zeros(2*model.mesh.nv,2*model.mesh.nv);
Zpp = zeros(model.mesh.np,model.mesh.np);
Zpt = zeros(model.mesh.np,4*model.mesh.nt);

if isunsteady == 1
    Gu = dt*Wi*H0M - dt*Wi*H1M - dt*2*BM';
    Gt = dt*C0M + dt*alpha*Wi*C1M + C2M;
else
    Gu = -Wi*H1M -2*BM';
    Gt = C0M + alpha*Wi*C1M;
end

% assemble full matrix
J = [Zvv -AM'   BM    ;
    AM    Zpp  Zpt  ;
    Gu   Zpt' Gt  ];

if isunsteady == 1
    Gu = -DV;
    Gt = -dt*alpha*Wi*C1V + Wi*C2V - dt*Wi*H0V + dt*Wi*H1V;
else
    Gu = -DV;
    Gt = -alpha*Wi*C1V + Wi*H1V;
end

% assemble vector
a = zeros(2*model.mesh.nv+model.mesh.np+4*model.mesh.nt,1);
a(1:2*model.mesh.nv) = Gu;
a(2*model.mesh.nv+model.mesh.np+1:2*model.mesh.nv+model.mesh.np+4*model.mesh.nt) = Gt;