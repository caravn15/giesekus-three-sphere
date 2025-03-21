function [I1,I2,I3,I4] = ElementalStiffnessH0Fixed(shape,dtN11dx,...
    dtN12dx,dtN21dx,dtN22dx,dtN11dy,dtN12dy,dtN21dy,...
        dtN22dy,uN,vN,detJ,w2D)

I1 = shape.t(:,:)*(dtN11dx.*uN.*(detJ.*w2D)) + ...
    shape.t(:,:)*(dtN11dy.*vN.*(detJ.*w2D));

I2 = shape.t(:,:)*(dtN12dx.*uN.*(detJ.*w2D)) + ...
    shape.t(:,:)*(dtN12dy.*vN.*(detJ.*w2D));

I3 = shape.t(:,:)*(dtN21dx.*uN.*(detJ.*w2D)) + ...
    shape.t(:,:)*(dtN21dy.*vN.*(detJ.*w2D));

I4 = shape.t(:,:)*(dtN22dx.*uN.*(detJ.*w2D)) + ...
    shape.t(:,:)*(dtN22dy.*vN.*(detJ.*w2D));

end