function [Jref,detJ] = CalculateRefJacobian(mesh)

Jref = cell(mesh.nelem,1);
detJ = cell(mesh.nelem,1);

for e=1:mesh.nelem
    
    % get element node coordinates
    xe = mesh.xv(mesh.vCon(e,:));
    ye = mesh.yv(mesh.vCon(e,:));
    
    % get limits in each direction
    x(e,:) = [min(xe),max(xe)];
    y(e,:) = [min(ye),max(ye)];

    % calculate element dimensions
    dx = x(e,2) - x(e,1);
    dy = y(e,2) - y(e,1);
    
    % calculate jacobian
    Jref{e} = [dx/2 0; 0 dy/2];
    
    % caculate determinant of jacobian
    detJ{e} = dx*dy/4;
    
end

end