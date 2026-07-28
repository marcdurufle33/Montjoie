function U = deriveX(X, Y, V)

xmin = min(min(X));
xmax = max(max(X));
nbx = size(X, 1);
nby = size(X, 2);
dy = (xmax-xmin)/(nby-1);
U = zeros(nbx, nby);
U(:, 3:nby-2) = (-V(:, 5:nby) + 8*V(:,4:nby-1) - 8*V(:, 2:nby-3) + V(:, 1:nby-4))/(12.0*dy);

U(:, 2) = (V(:, 3) - V(:, 1))/(2.0*dy);
U(:, nby-1) = (V(:, nby) - V(:, nby-2))/(2.0*dy);
U(:, 1) = (-3.0*V(:, 1) + 4.0*V(:, 2) - V(:, 3))/(2.0*dy);
U(:, nby) = (3.0*V(:, nby) - 4.0*V(:, nby-1) + V(:, nby-2))/(2.0*dy);

