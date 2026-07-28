function sor=plot2dinst(X, Y, M, cmin, cmax, xmin, xmax, ymin, ymax, num_figure)
% Fonction permettant d'afficher les instantannes.
% Necessite en entree:
%         1. matrices coordonnees: X et Y 
%         2. valeurs: 
%         3. les valeurs extremes a considerer (pour avoir un affichage a la
%         meme echelle de couleur sur les differents instantannes
%         4. on definit les axes: [xmin xmax]x[ymin ymax]
if(nargin==4)
  num_figure = cmin;
  cmin=min(min(M));
  cmax=max(max(M));
  xmin=min(min(X));
  ymin=min(min(Y));
  xmax=max(max(X));
  ymax=max(max(Y));
 
elseif (nargin==3)
  cmin=min(min(M));
  cmax=max(max(M));
  xmin=min(min(X));
  ymin=min(min(Y));
  xmax=max(max(X));
  ymax=max(max(Y));
  num_figure = 0;
elseif (nargin==5)
  xmin=min(min(X));
  ymin=min(min(Y));
  xmax=max(max(X));
  ymax=max(max(Y));
  num_figure = 0;
elseif (nargin==6)
  num_figure = xmin;
  xmin=min(min(X));
  ymin=min(min(Y));
  xmax=max(max(X));
  ymax=max(max(Y));
elseif (nargin==9)
  num_figure = 0;
end
format long
%hold on;
if (num_figure == 0)
  h = figure;
 % set(h,'paperposition', [0.25 2.5 5.675 4.95])
else
  %clf(num_figure);
end
%h = subplot('position',[.05 .55 .85 .4]);
%h = subplot('position',[.05 .05 .85 .4]);
hh = pcolor(X,Y,M);
shading interp;
caxis([cmin cmax])
colormap(jet(2000));
%colorbar
axis image
axis ([xmin xmax ymin ymax])
% drawnow
