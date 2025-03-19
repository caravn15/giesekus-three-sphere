function dshape = CalculateRefShapeDeriv(S)

% extract components
[s,t] = ExtractVectorComponents(S);

% velocity s derivatives
dshape.v(1,1,:) = 1/4*(2*s-1).*(t-1).*t;
dshape.v(1,2,:) = -s.*(t-1).*t;
dshape.v(1,3,:) = 1/4*(2*s+1).*(t-1).*t;
dshape.v(1,4,:) = -1/2*(2*s+1).*(t.^2-1);
dshape.v(1,5,:) = 1/4*(2*s+1).*t.*(t+1);
dshape.v(1,6,:) = -s.*t.*(t+1);
dshape.v(1,7,:) = 1/4*(2*s-1).*t.*(t+1);
dshape.v(1,8,:) = -1/2*(2*s-1).*(t.^2-1);
dshape.v(1,9,:) = 2*s.*(t.^2-1);

% velocity t derivatives
dshape.v(2,1,:) = 1/4*(s-1).*s.*(2*t-1);
dshape.v(2,2,:) = -1/2.*(s.^2-1).*(2*t-1);
dshape.v(2,3,:) = 1/4*s.*(s+1).*(2*t-1);
dshape.v(2,4,:) = -s.*(s+1).*t;
dshape.v(2,5,:) = 1/4*s.*(s+1).*(2*t+1);
dshape.v(2,6,:) = -1/2.*(s.^2-1).*(2*t+1);
dshape.v(2,7,:) = 1/4*(s-1).*s.*(2*t+1);
dshape.v(2,8,:) = -(s-1).*s.*t;
dshape.v(2,9,:) = 2*(s.^2-1).*t;

end