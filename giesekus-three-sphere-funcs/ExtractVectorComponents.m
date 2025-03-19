function [x1,x2] = ExtractVectorComponents(X)

N  = length(X)/2;
x1 = X(1:N);
x2 = X(N+1:2*N);

end 