function [I11x11,I11x12,I11x21,I11x22,I12x11,I12x12,I12x21,I12x22,...
        I21x11,I21x12,I21x21,I21x22,I22x11,I22x12,I22x21,I22x22] = ...
    ElementalStiffnessC1Mat(TpItQuad,shape,detJ,w2D)

[T11,T12,T21,T22] = ExtractTensorComponents(TpItQuad);

I11x11 = shape.t(:,:)*(shape.t(:,:)'.*(T11(:).*detJ.*w2D));

I12x11 = zeros(size(I11x11));

I21x11 = shape.t(:,:)*(shape.t(:,:)'.*(T21(:).*detJ.*w2D));

I22x11 = zeros(size(I11x11));

I11x12 = zeros(size(I11x11));

I12x12 = shape.t(:,:)*(shape.t(:,:)'.*(T11(:).*detJ.*w2D));

I21x12 = zeros(size(I11x11));

I22x12 = shape.t(:,:)*(shape.t(:,:)'.*(T21(:).*detJ.*w2D));

I11x21 = shape.t(:,:)*(shape.t(:,:)'.*(T12(:).*detJ.*w2D));

I12x21 = zeros(size(I11x11));

I21x21 = shape.t(:,:)*(shape.t(:,:)'.*(T22(:).*detJ.*w2D));

I22x21 = zeros(size(I11x11));

I11x22 = zeros(size(I11x11));

I12x22 = shape.t(:,:)*(shape.t(:,:)'.*(T12(:).*detJ.*w2D));

I21x22 = zeros(size(I11x11));

I22x22 = shape.t(:,:)*(shape.t(:,:)'.*(T22(:).*detJ.*w2D));

end