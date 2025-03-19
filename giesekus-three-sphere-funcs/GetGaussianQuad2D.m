function [x,w] = GetGaussianQuad2D(nquad)

% get one-dimensional gaussian quadrature rule on [-1,1]
[x1D,w1D]  = GetGaussianQuad1D(nquad);

% calculate two-dimensional GQ nodes
[x1,x2] = meshgrid(x1D,x1D);
x(:,1)  = x1(:);
x(:,2)  = x2(:);

% calculate two-dimensional GQ weights
[w1,w2] = meshgrid(w1D,w1D);
w       = w1(:).*w2(:);
w       = w(:);

end