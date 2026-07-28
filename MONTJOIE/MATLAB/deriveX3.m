function U = deriveX3(X, Y, Z, V)

xmin = min(X);
xmax = max(X);
nbx = length(X);
nby = length(Y);
nbz = length(Z);
dx = (xmax-xmin)/(nbx-1);
U = zeros(nbx, nby, nbz);
U(3:nbx-2, :, :) = 1.0/(12.0*dx)*(-V(5:nbx, :, :) + 8*V(4:nbx-1, :, :) - 8*V(2:nbx-3, :, :) + V(1:nbx-4, :, :));

U(2, :, :) = (V(3, :, :) - V(1, :, :))/(2.0*dx);
U(nbx-1, :, :) = (V(nbx, :, :) - V(nbx-2, :, :))/(2.0*dx);
U(1, :, :) = (-3.0*V(1, :, :) + 4.0*V(2, :, :) - V(3, :, :))/(2.0*dx);
U(nbx, :, :) = (3.0*V(nbx, :, :) - 4.0*V(nbx-1, :, :) + V(nbx-2, :, :))/(2.0*dx);

