function [x1,x2,x3,x4] = ExtractTensorComponents(X)

N  = length(X)/4;
x1 = X(1:N);
x2 = X(N+1:2*N);
x3 = X(2*N+1:3*N);
x4 = X(3*N+1:4*N);

end 