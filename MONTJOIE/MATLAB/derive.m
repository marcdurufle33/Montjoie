function U = laplacien(X, Y, V)

xmin = min(min(X));
xmax = max(max(X));
ymin = min(min(Y));
ymax = max(max(Y));
nbx = size(X, 1);
nby = size(Y, 2);
dx = (xmax-xmin)/(nbx-1)
dy = (ymax-ymin)/(nby-1)
U = zeros(nbx, nby);
U(2:nbx-1, :) = (-2*V(2:nbx-1, :) + V(3:nbx, :) + V(1:nbx-2, :))/(dx*dx);
U(1, :) = (2*V(1,:) - 5*V(2,:) + 4*V(3,:)-V(4,:))/(dx*dx);
U(nbx, :) = (2*V(nbx,:) - 5*V(nbx-1,:) + 4*V(nbx-2,:)-V(nbx-3,:))/(dx*dx);


U(:, 2:nby-1) = U(:, 2:nby-1) + (-2*V(:, 2:nby-1) + V(:, 3:nby) - V(:, 1:nby-2))/(dy*dy);
U(:, 1) = U(:, 1) + (2*V(:, 1) - 5*V(:, 2) + 4*V(:, 3)-V(:, 4))/(dy*dy);
U(:, nby) = U(:, nby) + (2*V(:, nby) - 5*V(:, nby-1) + 4*V(:, nby-2)-V(:, nby-3))/(dy*dy);
