function [I1x11,I1x12,I1x21,I1x22,I2x11,I2x12,I2x21,I2x22] = ...
    ElementalStiffnessC(TQuad,gdshape,shape,detJ,w2D)

[T11,T12,T21,T22] = ExtractTensorComponents(TQuad);

I1x11 = shape.t(:,:)*(squeeze(gdshape(1,:,:))'.*(T11(:).*detJ.*w2D)) + ...
    shape.t(:,:)*(squeeze(gdshape(1,:,:))'.*(T11(:).*detJ.*w2D));

I1x12 = shape.t(:,:)*(squeeze(gdshape(2,:,:))'.*(T11(:).*detJ.*w2D)) + ...
    shape.t(:,:)*(squeeze(gdshape(1,:,:))'.*(T12(:).*detJ.*w2D));

I1x21 = shape.t(:,:)*(squeeze(gdshape(1,:,:))'.*(T21(:).*detJ.*w2D)) + ...
    shape.t(:,:)*(squeeze(gdshape(2,:,:))'.*(T11(:).*detJ.*w2D));

I1x22 = shape.t(:,:)*(squeeze(gdshape(2,:,:))'.*(T21(:).*detJ.*w2D)) + ...
    shape.t(:,:)*(squeeze(gdshape(2,:,:))'.*(T12(:).*detJ.*w2D));

I2x11 = shape.t(:,:)*(squeeze(gdshape(1,:,:))'.*(T12(:).*detJ.*w2D)) + ...
    shape.t(:,:)*(squeeze(gdshape(1,:,:))'.*(T21(:).*detJ.*w2D));

I2x12 = shape.t(:,:)*(squeeze(gdshape(2,:,:))'.*(T12(:).*detJ.*w2D)) + ...
    shape.t(:,:)*(squeeze(gdshape(1,:,:))'.*(T22(:).*detJ.*w2D));

I2x21 = shape.t(:,:)*(squeeze(gdshape(1,:,:))'.*(T22(:).*detJ.*w2D)) + ...
    shape.t(:,:)*(squeeze(gdshape(2,:,:))'.*(T21(:).*detJ.*w2D));

I2x22 = shape.t(:,:)*(squeeze(gdshape(2,:,:))'.*(T22(:).*detJ.*w2D)) + ...
    shape.t(:,:)*(squeeze(gdshape(2,:,:))'.*(T22(:).*detJ.*w2D));

end

