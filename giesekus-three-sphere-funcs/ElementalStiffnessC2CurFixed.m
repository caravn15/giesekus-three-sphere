function [I1,I2,I3,I4] = ElementalStiffnessC2CurFixed(shape,tN11,tN12,...
    tN21,tN22,detJ,w2D)

I1 = shape.t(:,:)*(tN11.*(detJ.*w2D));

I2 = shape.t(:,:)*(tN12.*(detJ.*w2D));

I3 = shape.t(:,:)*(tN21.*(detJ.*w2D));

I4 = shape.t(:,:)*(tN22.*(detJ.*w2D));

end