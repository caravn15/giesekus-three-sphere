function U = GetFemVelocityAtPoint(mesh,U0,X)

% extract components
[x,y] = ExtractVectorComponents(X);

% get element index that point belongs to
e = GetFemElement(mesh,X);

% initilise
M = length(e);
U = zeros(2*M,1);
m = 1;

for ii=1:length(x)

    % get element node coordinates
    xe = mesh.xv(mesh.vCon(e(ii),:));
    ye = mesh.yv(mesh.vCon(e(ii),:));
    
    % get limits in each direction
    xlims(ii,:) = [min(xe),max(xe)];
    ylims(ii,:) = [min(ye),max(ye)];
    
    % transform point to reference element
    s = (2*x(ii)-(xlims(ii,1)+xlims(ii,2)))/(xlims(ii,2)-xlims(ii,1));
    t = (2*y(ii)-(ylims(ii,1)+ylims(ii,2)))/(ylims(ii,2)-ylims(ii,1));
    S = [s;t];

    % calculate value of shape functions at point
    shape = CalculateRefShape(S);

    % calculate velocities at nodes on element
    [Ux,Uy] = ExtractVectorComponents(U0);
    u       = Ux(mesh.vCon(e(ii),:));
    v       = Uy(mesh.vCon(e(ii),:));

    % interpolant to get fem velocity at point
    U(m)     = shape.v'*u;
    U(M+m)   = shape.v'*v;
    m = m+1;
end

end