function sor = film1D(V, dt)

sor = [];
N = size(V, 1);
x = linspace(-1, 1, N);
Vmin = min(min(V));
Vmax = max(max(V));
deltaV = Vmax - Vmin;
Vmin = Vmin - 0.1*deltaV;
Vmax = Vmax + 0.1*deltaV;
for i = 1:size(V,2)
  
  plot(x, V(:,i));
  axis([-1 1 Vmin Vmax]);
  drawnow;
  pause(dt)

end

