function dS = RegStokesletDerivVelocity(X,Y,ep)

[X1,X2] = ExtractVectorComponents(X);
[Y1,Y2] = ExtractVectorComponents(Y);

r1    = X1 - Y1';
r2    = X2 - Y2';
rsq   = r1.^2+r2.^2;
reps  = sqrt(rsq+ep^2);

dS{1,1,1} = (-2*(3*(r1.^2 - reps.^2/2 + ep^2)*ep.*sqrt(reps.^2) + ep^4 + (r1.^2 + reps.^2)*ep^2 + r1.^2.*reps.^2 - reps.^4/2).*r1./(3*(reps.^2).^(5/2)*ep + (reps.^2).^(3/2)*ep^3 + reps.^6 + 3*reps.^4*ep^2))/(4*pi);
dS{1,1,2} = (-2*(3*(r1.^2 + reps.^2/2 + ep^2)*ep.*sqrt(reps.^2) + ep^4 + (r1.^2 + 3*reps.^2)*ep^2 + r1.^2.*reps.^2 + reps.^4/2).*r2./(3*(reps.^2).^(5/2)*ep + (reps.^2).^(3/2)*ep^3 + reps.^6 + 3*reps.^4*ep^2))/(4*pi);
dS{1,2,1} = (-2*r2.*(3.*(r1.^2 - reps.^2/2)*ep.*sqrt(reps.^2) + (r1.^2 - reps.^2)*ep^2 + r1.^2.*reps.^2 - reps.^4/2)./(3*(reps.^2).^(5/2)*ep + (reps.^2).^(3/2)*ep^3 + reps.^6 + 3*reps.^4*ep^2))/(4*pi);
dS{1,2,2} = (2*r1.*(3*ep*(r1.^2 - reps.^2/2 + ep^2).*sqrt(reps.^2) + ep^4 + (r1.^2 + reps.^2)*ep^2 + r1.^2.*reps.^2 - reps.^4/2)./(3*(reps.^2).^(5/2)*ep + (reps.^2).^(3/2)*ep^3 + reps.^6 + 3*reps.^4*ep^2))/(4*pi);
dS{2,2,1} = (2*(3*(r1.^2 - (3*reps.^2)/2)*ep.*sqrt(reps.^2) + (r1.^2 - 3*reps.^2)*ep^2 + r1.^2.*reps.^2 - (3*reps.^4)/2).*r1./(3*(reps.^2).^(5/2)*ep + (reps.^2).^(3/2)*ep^3 + reps.^6 + 3*reps.^4*ep^2))/(4*pi);
dS{2,2,2} = (2*r2.*(3*ep*(r1.^2 - reps.^2/2).*sqrt(reps.^2) + (r1.^2 - reps.^2)*ep^2 + r1.^2.*reps.^2 - reps.^4/2)./(3*(reps.^2).^(5/2)*ep + (reps.^2).^(3/2)*ep^3 + reps.^6 + 3*reps.^4*ep^2))/(4*pi);
dS{2,1,1} = dS{1,2,1};
dS{2,1,2} = dS{1,2,2};

dS = gdivide(dS,4*pi);

end