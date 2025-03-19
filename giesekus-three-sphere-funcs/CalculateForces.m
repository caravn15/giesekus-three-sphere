function ft = CalculateForces(model,Xt,t)

% unpack model parameters
lr    = model.lr;
Amp   = model.Amp;
chi   = model.chi;
kappa = model.kappa;

% calculate spring length vectors
sl1 = Xt(:,2) - Xt(:,1);
sl2 = Xt(:,3) - Xt(:,2);

% calculate direction vectors
d1 = sl1/norm(sl1);
d2 = sl2/norm(sl2);

% calculate forces (hookean)
ftH(:,1) = (norm(sl1)-lr)*d1;
ftH(:,2) = -(norm(sl1)-lr)*d1 + (norm(sl2)-lr)*d2;
ftH(:,3) = -(norm(sl2)-lr)*d2;

% calculate forces (actuation)
ftA(:,1) = Amp*sin(2*pi*t)*d1;
ftA(:,2) = -Amp*sin(2*pi*t)*d1 + Amp*sin(2*pi*t-chi)*d2;
ftA(:,3) = -Amp*sin(2*pi*t-chi)*d2;

% calculate total forces
ft = kappa*(ftH+ftA);

end
