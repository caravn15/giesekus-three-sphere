function Xe = GetFemElement(mesh,X)

[x,y] = ExtractVectorComponents(X);

Xe = zeros(length(x),1);
for e=1:mesh.nelem

    % get element node coordinates
    xe = mesh.xv(mesh.vCon(e,:));
    ye = mesh.yv(mesh.vCon(e,:));

    % get limits in each direction
    xlims(e,:) = [min(xe),max(xe)];
    ylims(e,:) = [min(ye),max(ye)];

    x1 = (x >= xlims(e,1));
    x2 = (x <= xlims(e,2));
    y1 = (y >= ylims(e,1));
    y2 = (y <= ylims(e,2));

    Xt = x1.*x2.*y1.*y2;

    ind = find(Xt == 1);
    for ii=1:length(ind)
        if Xe(ind(ii)) == 0
            Xe(ind(ii)) = e;
        end
    end

    if(nnz(Xe) == length(Xe))
        break;
    end
end

end
