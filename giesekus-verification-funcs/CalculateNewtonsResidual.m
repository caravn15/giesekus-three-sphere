function [a,uAn] = CalculateNewtonsResidual(mesh,params)

%% calculate analytic solution

uAn = CalculateAnalyticSolution(params,mesh.lx(2),...
    mesh.yv(mesh.vbind{1}),'Giesekus');

%% construct full vector
a = AssembleResidualVector(mesh);

%% apply boundary conditions
a = ApplyBCsResidual(mesh,a,uAn);

end