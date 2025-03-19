function [x,w] = GetGaussianQuad1D(nQ)

a = -1;
b = 1;
nQ = nQ-1;
N1 = nQ+1; N2=nQ+2;
xu = linspace(-1,1,N1)';

% initial guess
y = cos((2*(0:nQ)'+1)*pi/(2*nQ+2))+(0.27/N1)*sin(pi*xu*nQ/N2);

% Legendre-Gauss Vandermonde matrix
L = zeros(N1,N2);

% derivative of LGVM
Lp = zeros(N1,N2);

% compute the zeros of the N+1 Legendre Polynomial
% using the recursion relation and the Newton-Raphson method
y0 = 2;

% Iterate until new points are uniformly within epsilon of old points
while max(abs(y-y0))>eps
    
    L(:,1)  = 1;
    Lp(:,1) = 0;
    
    L(:,2)  = y;
    Lp(:,2) = 1;
    
    for k=2:N1
        L(:,k+1) = ( (2*k-1)*y.*L(:,k)-(k-1)*L(:,k-1) )/k;
    end
 
    Lp = (N2)*( L(:,N1)-y.*L(:,N2) )./(1-y.^2);   
    
    y0 = y;
    y  = y0-L(:,N2)./Lp;
    
end

% linear map from[-1,1] to [a,b]
x = (a*(1-y)+b*(1+y))/2;  

% compute the weights
w = (b-a)./((1-y.^2).*Lp.^2)*(N2/N1)^2;

end