function a = ApplyBCsResidual(mesh,a,bc)

% separate blocks of residual vector corresponding to each component
au1 = a(1:mesh.nv);
au2 = a(mesh.nv+1:2*mesh.nv);
ap  = a(2*mesh.nv+1:2*mesh.nv+mesh.np);
at1 = a(2*mesh.nv+mesh.np+1:2*mesh.nv+mesh.np+mesh.nt);
at2 = a(2*mesh.nv+mesh.np+mesh.nt+1:2*mesh.nv+mesh.np+2*mesh.nt);
at3 = a(2*mesh.nv+mesh.np+2*mesh.nt+1:2*mesh.nv+mesh.np+3*mesh.nt);
at4 = a(2*mesh.nv+mesh.np+3*mesh.nt+1:2*mesh.nv+mesh.np+4*mesh.nt);

au1(bc.inds) = -bc.uN;
au2(bc.inds) = -bc.vN;

% reform residual vector
a = [au1; au2; ap; at1; at2; at3; at4];

end