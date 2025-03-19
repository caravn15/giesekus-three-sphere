function d2S = RegStokesletSecondDerivVelocity(X,Y,ep)

[X1,X2] = ExtractVectorComponents(X);
[Y1,Y2] = ExtractVectorComponents(Y);

r1    = X1 - Y1';
r2    = X2 - Y2';
rsq   = r1.^2+r2.^2;
reps  = sqrt(rsq+ep^2);

d2S{1,1,1,1} = (((24 .* r1 .^ 2 - 8 .* reps .^ 2) .* ep .^ 4 + (24 .* r1 .^ 4 ...
    - 16 .* r1 .^ 2 .* reps .^ 2 + reps .^ 4) .* ep .^ 2 + ...
    8 .* r1 .^ 4 .* reps .^ 2 - 8 .* r1 .^ 2 .* reps .^ 4 + reps .^ 6) ...
    .* sqrt((reps .^ 2)) + ((6 .* r1 .^ 2 - 2 .* reps .^ 2) .* ep .^ 5) ...
    + ((6 .* r1 .^ 4 + 26 .* r1 .^ 2 .* reps .^ 2 - 8 .* reps .^ 4) .* ep .^ 3) ...
    + ((32 .* r1 .^ 4 .* reps .^ 2 - 32 .* r1 .^ 2 .* reps .^ 4 + ...
    4 .* reps .^ 6) .* ep)) ./ ((reps .^ 2) .^ (0.9e1 ./ 0.2e1) + 6 ...
    .* (reps .^ 2) .^ (0.7e1 ./ 0.2e1) .* ep .^ 2 + ...
    (reps .^ 2) .^ (0.5e1 ./ 0.2e1) .* ep .^ 4 + 4 .* reps .^ 8 .* ep + ...
    4 .* reps .^ 6 .* ep .^ 3);

d2S{1,1,1,2} = (0.8e1 .* r2 .* ((3 .* ep .^ 4) + ((3 .* r1 .^ 2 + reps .^ 2) ...
    .* ep .^ 2) + (r1 .^ 2 .* reps .^ 2) - (reps .^ 4) ./ 0.4e1) .* r1 .* ...
    sqrt((reps .^ 2)) + 0.32e2 .* r2 .* ep .* r1 .* (0.3e1 ./ 0.16e2 .* ...
    (ep .^ 4) + (0.3e1 ./ 0.16e2 .* (r1 .^ 2) + (reps .^ 2)) .* (ep .^ 2) + ...
    (r1 .^ 2 .* reps .^ 2) - (reps .^ 4) ./ 0.4e1)) ./ ((reps .^ 2) .^ ...
    (0.9e1 ./ 0.2e1) + 6 .* (reps .^ 2) .^ (0.7e1 ./ 0.2e1) .* ep .^ 2 + ...
    (reps .^ 2) .^ (0.5e1 ./ 0.2e1) .* ep .^ 4 + 4 .* reps .^ 8 .* ep + ...
    4 .* reps .^ 6 .* ep .^ 3);

d2S{1,1,2,1} = d2S{1,1,1,2};

d2S{1,1,2,2} = ((-24 .* ep .^ 6 + (-48 .* r1 .^ 2 - 8 .* reps .^ 2) .* ...
    ep .^ 4 + (-24 .* r1 .^ 4 - 16 .* r1 .^ 2 .* reps .^ 2 + 13 .* reps .^ 4) ...
    .* ep .^ 2 - 8 .* r1 .^ 4 .* reps .^ 2 + 4 .* r1 .^ 2 .* reps .^ 4 + reps .^ ...
    6) .* sqrt((reps .^ 2)) - (6 .* ep .^ 7) + ((-12 .* r1 .^ 2 - ...
    32 .* reps .^ 2) .* ep .^ 5) + ((-6 .* r1 .^ 4 - 64 .* r1 .^ 2 .* reps .^ ...
    2 + 16 .* reps .^ 4) .* ep .^ 3) + ((-32 .* r1 .^ 4 .* reps .^ 2 + ...
    16 .* r1 .^ 2 .* reps .^ 4 + 4 .* reps .^ 6) .* ep)) ./ ((reps .^ 2) .^ ...
    (0.9e1 ./ 0.2e1) + 6 .* (reps .^ 2) .^ (0.7e1 ./ 0.2e1) .* ep .^ 2 + ...
    (reps .^ 2) .^ (0.5e1 ./ 0.2e1) .* ep .^ 4 + 4 .* reps .^ 8 .* ep + 4 .* ...
    reps .^ 6 .* ep .^ 3);

d2S{1,2,1,1} = 0.2e1 ./ ((reps .^ 2) .^ (0.9e1 ./ 0.2e1) + 6 .* ...
    (reps .^ 2) .^ (0.7e1 ./ 0.2e1) .* ep .^ 2 + (reps .^ 2) .^ ...
    (0.5e1 ./ 0.2e1) .* ep .^ 4 + 4 .* reps .^ 8 .* ep + 4 .* reps .^ ...
    6 .* ep .^ 3) .* r1 .* r2 .* (0.4e1 .* sqrt((reps .^ 2)) .* r1 .^ 2 .* ...
    (reps .^ 2) + 0.12e2 .* sqrt((reps .^ 2)) .* r1 .^ 2 .* (ep .^ 2) - ...
    0.3e1 .* sqrt((reps .^ 2)) .* (reps .^ 4) - 0.12e2 .* sqrt((reps .^ 2)) ...
    .* (reps .^ 2) .* (ep .^ 2) + 0.16e2 .* r1 .^ 2 .* (reps .^ 2) .* ep + ...
    0.3e1 .* (ep .^ 3) .* r1 .^ 2 - (12 .* reps .^ 4 .* ep) - ...
    (3 .* reps .^ 2 .* ep .^ 3));

d2S{1,2,1,2} = (((-24 .* r1 .^ 2 + 8 .* reps .^ 2) .* ep .^ 4 + (-24 .* r1 .^ ...
    4 + 16 .* r1 .^ 2 .* reps .^ 2 - reps .^ 4) .* ep .^ 2 - 8 .* r1 .^ 4 .* ...
    reps .^ 2 + 8 .* r1 .^ 2 .* reps .^ 4 - reps .^ 6) .* sqrt((reps .^ 2)) + ...
    ((-6 .* r1 .^ 2 + 2 .* reps .^ 2) .* ep .^ 5) + ((-6 .* r1 .^ 4 - 26 .* r1 ...
    .^ 2 .* reps .^ 2 + 8 .* reps .^ 4) .* ep .^ 3) + ((-32 .* r1 .^ 4 .* ...
    reps .^ 2 + 32 .* r1 .^ 2 .* reps .^ 4 - 4 .* reps .^ 6) .* ep)) ./ ...
    ((reps .^ 2) .^ (0.9e1 ./ 0.2e1) + 6 .* (reps .^ 2) .^ (0.7e1 ./ 0.2e1) ...
    .* ep .^ 2 + (reps .^ 2) .^ (0.5e1 ./ 0.2e1) .* ep .^ 4 + 4 .* reps .^ 8 .* ...
    ep + 4 .* reps .^ 6 .* ep .^ 3);

d2S{1,2,2,1} = d2S{1,2,1,2};

d2S{1,2,2,2} = (-0.8e1 .* r2 .* ((3 .* ep .^ 4) + ((3 .* r1 .^ 2 + reps .^ ...
    2) .* ep .^ 2) + (r1 .^ 2 .* reps .^ 2) - (reps .^ 4) ./ 0.4e1) .* r1 .* ...
    sqrt((reps .^ 2)) - 0.32e2 .* r2 .* ep .* r1 .* (0.3e1 ./ 0.16e2 .* ...
    (ep .^ 4) + (0.3e1 ./ 0.16e2 .* (r1 .^ 2) + (reps .^ 2)) .* (ep .^ 2) + ...
    (r1 .^ 2 .* reps .^ 2) - (reps .^ 4) ./ 0.4e1)) ./ ((reps .^ 2) .^ ...
    (0.9e1 ./ 0.2e1) + 6 .* (reps .^ 2) .^ (0.7e1 ./ 0.2e1) .* ep .^ 2 + ...
    (reps .^ 2) .^ (0.5e1 ./ 0.2e1) .* ep .^ 4 + 4 .* reps .^ 8 .* ep + 4 .* ...
    reps .^ 6 .* ep .^ 3);

d2S{2,1,1,1} = d2S{1,2,1,1};

d2S{2,1,1,2} = d2S{1,2,1,2};

d2S{2,1,2,1} = d2S{1,2,2,1};

d2S{2,1,2,2} = d2S{1,2,2,2};

d2S{2,2,1,1} = (((-24 .* r1 .^ 4 + 48 .* r1 .^ 2 .* reps .^ 2 - 15 .* reps ...
    .^ 4) .* ep .^ 2 - 8 .* r1 .^ 4 .* reps .^ 2 + 12 .* r1 .^ 2 .* reps .^ 4 - ...
    3 .* reps .^ 6) .* sqrt((reps .^ 2)) - 0.32e2 .* ep .* (0.3e1 ./ 0.16e2 ...
    .* ((r1 - reps) .^ 2) .* ((r1 + reps) .^ 2) .* (ep .^ 2) + (reps .^ 2) .* ...
    ((r1 .^ 4) - 0.3e1 ./ 0.2e1 .* (r1 .^ 2) .* (reps .^ 2) + 0.3e1 ./ ...
    0.8e1 .* (reps .^ 4)))) ./ ((reps .^ 2) .^ (0.9e1 ./ 0.2e1) + 6 .* ...
    (reps .^ 2) .^ (0.7e1 ./ 0.2e1) .* ep .^ 2 + (reps .^ 2) .^ (0.5e1 ./ ...
    0.2e1) .* ep .^ 4 + 4 .* reps .^ 8 .* ep + 4 .* reps .^ 6 .* ep .^ 3);

d2S{2,2,1,2} = -0.2e1 ./ ((reps .^ 2) .^ (0.9e1 ./ 0.2e1) + 6 .* (reps .^ ...
    2) .^ (0.7e1 ./ 0.2e1) .* ep .^ 2 + (reps .^ 2) .^ (0.5e1 ./ 0.2e1) .* ep ...
    .^ 4 + 4 .* reps .^ 8 .* ep + 4 .* reps .^ 6 .* ep .^ 3) .* r1 .* r2 .* ...
    (0.4e1 .* sqrt((reps .^ 2)) .* r1 .^ 2 .* (reps .^ 2) + 0.12e2 .* ...
    sqrt((reps .^ 2)) .* r1 .^ 2 .* (ep .^ 2) - 0.3e1 .* sqrt((reps .^ 2)) ...
    .* (reps .^ 4) - 0.12e2 .* sqrt((reps .^ 2)) .* (reps .^ 2) .* (ep .^ 2) + ...
    0.16e2 .* r1 .^ 2 .* (reps .^ 2) .* ep + 0.3e1 .* (ep .^ 3) .* r1 .^ 2 - ...
    (12 .* reps .^ 4 .* ep) - (3 .* reps .^ 2 .* ep .^ 3));

d2S{2,2,2,1} = d2S{2,2,1,2};

d2S{2,2,2,2} = (((24 .* r1 .^ 2 - 8 .* reps .^ 2) .* ep .^ 4 + (24 .* r1 ...
    .^ 4 - 16 .* r1 .^ 2 .* reps .^ 2 + reps .^ 4) .* ep .^ 2 + 8 .* r1 .^ 4 .* ...
    reps .^ 2 - 8 .* r1 .^ 2 .* reps .^ 4 + reps .^ 6) .* sqrt((reps .^ 2)) + ...
    ((6 .* r1 .^ 2 - 2 .* reps .^ 2) .* ep .^ 5) + ((6 .* r1 .^ 4 + 26 .* r1 .^ ...
    2 .* reps .^ 2 - 8 .* reps .^ 4) .* ep .^ 3) + ((32 .* r1 .^ 4 .* reps .^ ...
    2 - 32 .* r1 .^ 2 .* reps .^ 4 + 4 .* reps .^ 6) .* ep)) ./ ((reps .^ 2) .^ ...
    (0.9e1 ./ 0.2e1) + 6 .* (reps .^ 2) .^ (0.7e1 ./ 0.2e1) .* ep .^ 2 + ...
    (reps .^ 2) .^ (0.5e1 ./ 0.2e1) .* ep .^ 4 + 4 .* reps .^ 8 .* ep + 4 .* ...
    reps .^ 6 .* ep .^ 3);

d2S = gdivide(d2S,4*pi);

end