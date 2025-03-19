function [I1,I2,I3,I4] = ElementalStiffnessC1(shape,tN11,tN12,...
    tN21,tN22,TpItQuad,detJ,w2D)

[T11,T12,T21,T22] = ExtractTensorComponents(TpItQuad);

I1 = shape.t(:,:)*(tN11.*T11.*(detJ.*w2D)) + ...
    + shape.t(:,:)*(tN21.*T12.*(detJ.*w2D));

I2 = shape.t(:,:)*(tN12.*T11.*(detJ.*w2D)) + ...
    shape.t(:,:)*(tN22.*T12.*(detJ.*w2D));

I3 = shape.t(:,:)*(tN11.*T21.*(detJ.*w2D)) + ...
    shape.t(:,:)*(tN21.*T22.*(detJ.*w2D));

I4 = shape.t(:,:)*(tN12.*T21.*(detJ.*w2D)) + ...
    shape.t(:,:)*(tN22.*T22.*(detJ.*w2D));

end