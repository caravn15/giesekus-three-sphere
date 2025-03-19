function shape = CalculateRefShape(S)

% extract components
[s,t] = ExtractVectorComponents(S);

% velocity shape functions
shape.v(1,:) = 1/4*s.*(s-1).*t.*(t-1);
shape.v(2,:) = -1/2.*(s+1).*(s-1).*t.*(t-1);
shape.v(3,:) = 1/4*(s+1).*s.*t.*(t-1);
shape.v(4,:) = -1/2*(s+1).*s.*(t+1).*(t-1);
shape.v(5,:) = 1/4*(s+1).*s.*(t+1).*t;
shape.v(6,:) = -1/2*(s+1).*(s-1).*(t+1).*t;
shape.v(7,:) = 1/4*s.*(s-1).*(t+1).*t;
shape.v(8,:) = -1/2*s.*(s-1).*(t+1).*(t-1);
shape.v(9,:) = (s+1).*(s-1).*(t+1).*(t-1);

% pressure shape functions
shape.p(1,:) = 1/4*(1-s).*(1-t);
shape.p(2,:) = 1/4*(1+s).*(1-t);
shape.p(3,:) = 1/4*(1+s).*(1+t);
shape.p(4,:) = 1/4*(1-s).*(1+t);

% stress shape functions
shape.t(1,:) = 1/4*s.*(s-1).*t.*(t-1);
shape.t(2,:) = -1/2.*(s+1).*(s-1).*t.*(t-1);
shape.t(3,:) = 1/4*(s+1).*s.*t.*(t-1);
shape.t(4,:) = -1/2*(s+1).*s.*(t+1).*(t-1);
shape.t(5,:) = 1/4*(s+1).*s.*(t+1).*t;
shape.t(6,:) = -1/2*(s+1).*(s-1).*(t+1).*t;
shape.t(7,:) = 1/4*s.*(s-1).*(t+1).*t;
shape.t(8,:) = -1/2*s.*(s-1).*(t+1).*(t-1);
shape.t(9,:) = (s+1).*(s-1).*(t+1).*(t-1);

end