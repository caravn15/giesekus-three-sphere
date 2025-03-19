function TQuad = GetStressAtQuadPoints(mesh,shape,taupr)

% extract stress components
[T11,T12,T21,T22] = ExtractTensorComponents(taupr);

for e=1:mesh.nelem
    
    % get stress components on element
    T11e = T11(mesh.tCon(e,:));
    T12e = T12(mesh.tCon(e,:));
    T21e = T21(mesh.tCon(e,:));
    T22e = T22(mesh.tCon(e,:));
    
    % get stress at quad points
    T11Quad = T11e'*shape.t(:,:);
    T12Quad = T12e'*shape.t(:,:);
    T21Quad = T21e'*shape.t(:,:);
    T22Quad = T22e'*shape.t(:,:);
    
    % sort output
    TQuad{e} = [T11Quad(:); T12Quad(:); T21Quad(:); T22Quad(:)];
    
end

end