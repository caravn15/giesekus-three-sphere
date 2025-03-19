function [I1x11,I1x12,I1x21,I1x22,I2x11,I2x12,I2x21,I2x22] = ...
    ElementalStiffnessB(gdshape,shape,detJ,w2D)

I1x11 = squeeze(gdshape(1,:,:))*(shape.t(:,:)'.*(detJ*w2D));
I1x12 = 1/2*squeeze(gdshape(2,:,:))*(shape.t(:,:)'.*(detJ*w2D));
I1x21 = 1/2*squeeze(gdshape(2,:,:))*(shape.t(:,:)'.*(detJ*w2D));
I1x22 = zeros(size(I1x11));

I2x11 = zeros(size(I1x11));
I2x12 = 1/2*squeeze(gdshape(1,:,:))*(shape.t(:,:)'.*(detJ*w2D));
I2x21 = 1/2*squeeze(gdshape(1,:,:))*(shape.t(:,:)'.*(detJ*w2D));
I2x22 = squeeze(gdshape(2,:,:))*(shape.t(:,:)'.*(detJ*w2D));

end