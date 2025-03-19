function [I1,I2] = ElementalStiffnessA(gdshape,shape,detJ,w2D)

I1 = shape.p(:,:)*(squeeze(gdshape(1,:,:))'.*(detJ*w2D));
I2 = shape.p(:,:)*(squeeze(gdshape(2,:,:))'.*(detJ*w2D));

end