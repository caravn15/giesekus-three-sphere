function [I1,I2,I3,I4] = ElementalStiffnessH0(shape,...
        dTpItQuad,uN,vN,detJ,w2D)

% get tensor components
[dTpIt11dx,dTpIt12dx,dTpIt21dx,dTpIt22dx] = ...
    ExtractTensorComponents(dTpItQuad.x);

[dTpIt11dy,dTpIt12dy,dTpIt21dy,dTpIt22dy] = ...
    ExtractTensorComponents(dTpItQuad.y);

I1 = shape.t(:,:)*(dTpIt11dx.*uN.*(detJ.*w2D)) + ...
    shape.t(:,:)*(dTpIt11dy.*vN.*(detJ.*w2D));

I2 = shape.t(:,:)*(dTpIt12dx.*uN.*(detJ.*w2D)) + ...
    shape.t(:,:)*(dTpIt12dy.*vN.*(detJ.*w2D));

I3 = shape.t(:,:)*(dTpIt21dx.*uN.*(detJ.*w2D)) + ...
    shape.t(:,:)*(dTpIt21dy.*vN.*(detJ.*w2D));

I4 = shape.t(:,:)*(dTpIt22dx.*uN.*(detJ.*w2D)) + ...
    shape.t(:,:)*(dTpIt22dy.*vN.*(detJ.*w2D));

end