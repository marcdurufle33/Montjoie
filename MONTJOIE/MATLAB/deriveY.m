function U = deriveY(X, Y, V)

xmin = min(min(Y));
xmax = max(max(Y));
nbx = size(X, 1);
nby = size(X, 2);
dx = (xmax-xmin)/(nby-1);
U = zeros(nbx, nby);
U(3:nbx-2,:) = (-V(5:nbx,:) + 8*V(4:nbx-1,:) - 8*V(2:nbx-3,:) + V(1:nbx-4,:))/(12.0*dx);
U(2, :) = (V(3, :) - V(1, :))/(2.0*dx);
U(nbx-1, :) = (V(nbx, :) - V(nbx-2, :))/(2.0*dx);
U(1, :) = (-3.0*V(1,:) + 4.0*V(2, :) - V(3, :))/(2.0*dx);
U(nbx, :) = (3.0*V(nbx, :) - 4.0*V(nbx-1, :) + V(nbx-2,:))/(2.0*dx);

