function [I11x1,I12x1,I21x1,I22x1,I11x2,I12x2,I21x2,I22x2] = ...
    ElementalStiffnessH1FixedMat(gdshape,shape,tN11,tN12,tN21,tN22,...
    detJ,w2D)

I11x1 = shape.t(:,:)*(squeeze(gdshape(1,:,:))'.*(tN11(:).*detJ.*w2D)) + ...
    shape.t(:,:)*(squeeze(gdshape(2,:,:))'.*(tN12(:).*detJ.*w2D)) + ...
    shape.t(:,:)*(squeeze(gdshape(1,:,:))'.*(tN11(:).*detJ.*w2D)) + ...
    shape.t(:,:)*(squeeze(gdshape(2,:,:))'.*(tN21(:).*detJ.*w2D));

I12x1 = shape.t(:,:)*(squeeze(gdshape(1,:,:))'.*(tN12(:).*detJ.*w2D)) + ...
    shape.t(:,:)*(squeeze(gdshape(2,:,:))'.*(tN22(:).*detJ.*w2D));

I21x1 = shape.t(:,:)*(squeeze(gdshape(1,:,:))'.*(tN21(:).*detJ.*w2D)) + ...
    shape.t(:,:)*(squeeze(gdshape(2,:,:))'.*(tN22(:).*detJ.*w2D));

I22x1 = zeros(size(I11x1));

I11x2 = zeros(size(I11x1));

I12x2 = shape.t(:,:)*(squeeze(gdshape(1,:,:))'.*(tN11(:).*detJ.*w2D)) + ...
    shape.t(:,:)*(squeeze(gdshape(2,:,:))'.*(tN12(:).*detJ.*w2D));

I21x2 = shape.t(:,:)*(squeeze(gdshape(1,:,:))'.*(tN11(:).*detJ.*w2D)) + ...
    shape.t(:,:)*(squeeze(gdshape(2,:,:))'.*(tN21(:).*detJ.*w2D));

I22x2 = shape.t(:,:)*(squeeze(gdshape(1,:,:))'.*(tN21(:).*detJ.*w2D)) + ...
    shape.t(:,:)*(squeeze(gdshape(2,:,:))'.*(tN22(:).*detJ.*w2D)) + ...
    shape.t(:,:)*(squeeze(gdshape(1,:,:))'.*(tN12(:).*detJ.*w2D)) + ...
    shape.t(:,:)*(squeeze(gdshape(2,:,:))'.*(tN22(:).*detJ.*w2D));

end

