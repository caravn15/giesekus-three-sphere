function [I11x1,I12x1,I21x1,I22x1,I11x2,I12x2,I21x2,I22x2] = ...
    ElementalStiffnessH0Mat(shape,dTpItQuad,...
    detJ,w2D)

% get tensor components
[dTpIt11dx,dTpIt12dx,dTpIt21dx,dTpIt22dx] = ...
    ExtractTensorComponents(dTpItQuad.x);

[dTpIt11dy,dTpIt12dy,dTpIt21dy,dTpIt22dy] = ...
    ExtractTensorComponents(dTpItQuad.y);

I11x1 = shape.t(:,:)*(shape.v(:,:)'.*(dTpIt11dx(:).*detJ.*w2D));

I12x1 = shape.t(:,:)*(shape.v(:,:)'.*(dTpIt12dx(:).*detJ.*w2D));

I21x1 = shape.t(:,:)*(shape.v(:,:)'.*(dTpIt21dx(:).*detJ.*w2D));

I22x1 = shape.t(:,:)*(shape.v(:,:)'.*(dTpIt22dx(:).*detJ.*w2D));

I11x2 = shape.t(:,:)*(shape.v(:,:)'.*(dTpIt11dy(:).*detJ.*w2D));

I12x2 = shape.t(:,:)*(shape.v(:,:)'.*(dTpIt12dy(:).*detJ.*w2D));

I21x2 = shape.t(:,:)*(shape.v(:,:)'.*(dTpIt21dy(:).*detJ.*w2D));

I22x2 = shape.t(:,:)*(shape.v(:,:)'.*(dTpIt22dy(:).*detJ.*w2D));

end

