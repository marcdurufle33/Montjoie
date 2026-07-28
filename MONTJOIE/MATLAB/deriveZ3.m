function U = deriveZ3(X, Y, Z, V)

zmin = min(Z);
zmax = max(Z);
nbx = length(X);
nby = length(Y);
nbz = length(Z);
dz = (zmax-zmin)/(nbz-1);
U = zeros(nbx, nby, nbz);
U(:, :, 3:nbz-2) = 1.0/(12.0*dz)*(-V(:, :, 5:nbz) + 8*V(:, :, 4:nbz-1) - 8*V(:, :, 2:nbz-3) + V(:, :, 1:nbz-4));

U(:, :, 2) = (V(:, :, 3) - V(:, :, 1))/(2.0*dz);
U(:, :, nbz-1) = (V(:, :, nbz) - V(:, :, nbz-2))/(2.0*dz);
U(:, :, 1) = (-3.0*V(:, :, 1) + 4.0*V(:, :, 2) - V(:, :, 3))/(2.0*dz);
U(:, :, nbz) = (3.0*V(:, :, nbz) - 4.0*V(:, :, nbz-1) + V(:, :, nbz-2))/(2.0*dz);

