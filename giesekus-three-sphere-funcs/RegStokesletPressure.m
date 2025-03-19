function P = RegStokesletPressure(X,Y,ep)

[X1,X2] = ExtractVectorComponents(X);
[Y1,Y2] = ExtractVectorComponents(Y);

r1   = X1 - Y1';
r2   = X2 - Y2';
rsq  = r1.^2+r2.^2;
reps = sqrt(rsq+ep^2);

A  = ((rsq+2*ep^2+ep*reps)./((reps+ep).*(reps.^2).^(3/2)));
P1 = r1.*A;
P2 = r2.*A;

P = 1/(2*pi)*[P1 P2];

end