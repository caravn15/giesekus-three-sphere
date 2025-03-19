function J = ApplyBCsJacobian(mesh,J)

% separate blocks of Jacobian corresponding to each component
Ju1 = J(1:mesh.nv,:);
Ju2 = J(mesh.nv+1:2*mesh.nv,:);
Jp  = J(2*mesh.nv+1:2*mesh.nv+mesh.np,:);
Jt1 = J(2*mesh.nv+mesh.np+1:2*mesh.nv+mesh.np+mesh.nt,:);
Jt2 = J(2*mesh.nv+mesh.np+mesh.nt+1:2*mesh.nv+mesh.np+2*mesh.nt,:);
Jt3 = J(2*mesh.nv+mesh.np+2*mesh.nt+1:2*mesh.nv+mesh.np+3*mesh.nt,:);
Jt4 = J(2*mesh.nv+mesh.np+3*mesh.nt+1:2*mesh.nv+mesh.np+4*mesh.nt,:);

% get velocity boundary indices from mesh structure
binds = [mesh.vbind{1}'; mesh.vbind{2}'; mesh.vbind{4}'];

% set velocity boundary rows to zero
Ju1(binds,:) = 0;
Ju2(binds,:) = 0;

% get matrix indicies of velocity boundary nodes
ind1  = sub2ind(size(Ju1),binds,binds);
ind2  = sub2ind(size(Ju1),binds,binds+mesh.nv);

% update to prescribe velocity at boundary nodes
Ju1(ind1) = 1;
Ju2(ind2) = 1;

% reform Jacobian matrix
J = [Ju1; Ju2; Jp; Jt1; Jt2; Jt3; Jt4];

end