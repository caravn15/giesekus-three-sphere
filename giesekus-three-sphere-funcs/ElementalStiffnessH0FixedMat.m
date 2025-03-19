function [I11x1,I12x1,I21x1,I22x1,I11x2,I12x2,I21x2,I22x2] = ...
    ElementalStiffnessH0FixedMat(shape,dtN11dx,...
    dtN12dx,dtN21dx,dtN22dx,dtN11dy,dtN12dy,dtN21dy,...
        dtN22dy,detJ,w2D)

I11x1 = shape.t(:,:)*(shape.v(:,:)'.*(dtN11dx(:).*detJ.*w2D));

I12x1 = shape.t(:,:)*(shape.v(:,:)'.*(dtN12dx(:).*detJ.*w2D));

I21x1 = shape.t(:,:)*(shape.v(:,:)'.*(dtN21dx(:).*detJ.*w2D));

I22x1 = shape.t(:,:)*(shape.v(:,:)'.*(dtN22dx(:).*detJ.*w2D));

I11x2 = shape.t(:,:)*(shape.v(:,:)'.*(dtN11dy(:).*detJ.*w2D));

I12x2 = shape.t(:,:)*(shape.v(:,:)'.*(dtN12dy(:).*detJ.*w2D));

I21x2 = shape.t(:,:)*(shape.v(:,:)'.*(dtN21dy(:).*detJ.*w2D));

I22x2 = shape.t(:,:)*(shape.v(:,:)'.*(dtN22dy(:).*detJ.*w2D));

end