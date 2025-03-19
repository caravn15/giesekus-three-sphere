function S = RegStokesletVelocity(X,Y,ep)

[X1,X2] = ExtractVectorComponents(X);
[Y1,Y2] = ExtractVectorComponents(Y);

r1   = X1 - Y1';
r2   = X2 - Y2';
rsq  = r1.^2+r2.^2;
reps = sqrt(rsq+ep^2);

A = -log(reps+ep)+(ep*(reps+2*ep))./((reps+ep).*reps);
B = (reps+2*ep)./((reps+ep).^2.*reps); 

S11 = A + r1.*r1.*B;
S12 = r1.*r2.*B;
S21 = r2.*r1.*B;
S22 = A+r2.*r2.*B;

S = 1/(4*pi)*[S11 S12; S21 S22];

end