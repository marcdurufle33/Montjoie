function x = plot_triangle(pos, order, taille, couleur)

x0 = pos(1);
y0 = pos(2);
dx = taille;
dy = taille*order;
t = line([x0, x0+dx], [y0 y0]);
set(t, 'Color', couleur);
t = line([x0+dx x0], [y0 y0-dy]);
set(t, 'Color', couleur);
t = line([x0, x0], [y0-dy y0]);
set(t, 'Color', couleur);

text(x0+0.4*dx, y0+0.1*dy, '1');
text(x0-0.1*dx, y0-0.6*dy, num2str(order));
