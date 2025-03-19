function [x,y] = TransformFromGlobal(y0t,x2D,model)

% unpack required parameters
ep   = model.ep;
xlim = model.mesh.xlim;
ylim = model.mesh.ylim;

% calculate m parameter in sinh transform
for jj=1:2
    if abs(y0t(jj))<1
        alpha = ep;
    else
        %         alpha = min(ep,abs(y0t(jj))-1);
        alpha = ep;
    end

    m(jj) = floor(1 + log(alpha^(-1)));
end


x = zeros(length(x2D),model.mesh.nelem);
y = zeros(length(x2D),model.mesh.nelem);
for e=1:model.mesh.nelem

    % calculate transformed integral limits
    xlimt = sinh(1/m(1)*asinh((xlim(e,:)-y0t(1))/ep));
    ylimt = sinh(1/m(2)*asinh((ylim(e,:)-y0t(2))/ep));

    % calculate transformed quadrature points
    X = (xlimt(2)-xlimt(1))/2*x2D(:,1)+(xlimt(2)+xlimt(1))/2;
    Y = (ylimt(2)-ylimt(1))/2*x2D(:,2)+(ylimt(2)+ylimt(1))/2;

    x(:,e) = y0t(1)+ep*sinh(m(1)*asinh(X));
    y(:,e) = y0t(2)+ep*sinh(m(2)*asinh(Y));
end

end