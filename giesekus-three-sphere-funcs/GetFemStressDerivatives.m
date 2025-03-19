function dT = GetFemStressDerivatives(Tp,mesh,gdshape,varargin)

% extract stress components
[T11,T12,T21,T22] = ExtractTensorComponents(Tp);

if length(varargin) == 1

    e = varargin{1};

    % get stress components on element
    T11e = T11(mesh.tCon(e,:));
    T12e = T12(mesh.tCon(e,:));
    T21e = T21(mesh.tCon(e,:));
    T22e = T22(mesh.tCon(e,:));

    % get stress x derivatives
    dt11dxGQ = T11e'*squeeze(gdshape(1,:,:));
    dt12dxGQ = T12e'*squeeze(gdshape(1,:,:));
    dt21dxGQ = T21e'*squeeze(gdshape(1,:,:));
    dt22dxGQ = T22e'*squeeze(gdshape(1,:,:));

    % get stress y derivatives
    dt11dyGQ = T11e'*squeeze(gdshape(2,:,:));
    dt12dyGQ = T12e'*squeeze(gdshape(2,:,:));
    dt21dyGQ = T21e'*squeeze(gdshape(2,:,:));
    dt22dyGQ = T22e'*squeeze(gdshape(2,:,:));

    % sort output
    dT.x = [dt11dxGQ(:); dt12dxGQ(:); dt21dxGQ(:); dt22dxGQ(:)];
    dT.y = [dt11dyGQ(:); dt12dyGQ(:); dt21dyGQ(:); dt22dyGQ(:)];

else

    for e=1:mesh.nelem

        % get stress components on element
        T11e = T11(mesh.tCon(e,:));
        T12e = T12(mesh.tCon(e,:));
        T21e = T21(mesh.tCon(e,:));
        T22e = T22(mesh.tCon(e,:));

        % get stress x derivatives
        dt11dxGQ = T11e'*squeeze(gdshape{e}(1,:,:));
        dt12dxGQ = T12e'*squeeze(gdshape{e}(1,:,:));
        dt21dxGQ = T21e'*squeeze(gdshape{e}(1,:,:));
        dt22dxGQ = T22e'*squeeze(gdshape{e}(1,:,:));

        % get stress y derivatives
        dt11dyGQ = T11e'*squeeze(gdshape{e}(2,:,:));
        dt12dyGQ = T12e'*squeeze(gdshape{e}(2,:,:));
        dt21dyGQ = T21e'*squeeze(gdshape{e}(2,:,:));
        dt22dyGQ = T22e'*squeeze(gdshape{e}(2,:,:));

        % sort output
        dT{e}.x = [dt11dxGQ(:); dt12dxGQ(:); dt21dxGQ(:); dt22dxGQ(:)];
        dT{e}.y = [dt11dyGQ(:); dt12dyGQ(:); dt21dyGQ(:); dt22dyGQ(:)];

    end
end
end