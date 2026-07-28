function U = deriveY3(X, Y, Z, V)

ymin = min(Y);
ymax = max(Y);
nbx = length(X);
nby = length(Y);
nbz = length(Z);
dy = (ymax-ymin)/(nby-1);
U = zeros(nbx, nby, nbz);
U(:, 3:nby-2, :) = 1.0/(12.0*dy)*(-V(:, 5:nby, :) + 8*V(:, 4:nby-1, :) - 8*V(:, 2:nby-3, :) + V(:, 1:nby-4, :));

U(:, 2, :) = (V(:, 3, :) - V(:, 1, :))/(2.0*dy);
U(:, nby-1, :) = (V(:, nby, :) - V(:, nby-2, :))/(2.0*dy);
U(:, 1, :) = (-3.0*V(:, 1, :) + 4.0*V(:, 2, :) - V(:, 3, :))/(2.0*dy);
U(:, nby, :) = (3.0*V(:, nby, :) - 4.0*V(:, nby-1, :) + V(:, nby-2, :))/(2.0*dy);

