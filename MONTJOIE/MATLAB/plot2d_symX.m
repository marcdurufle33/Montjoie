function sor=plot2d_symY(X,Y,M,cmin,cmax,xmin,xmax,ymin,ymax)
% Fonction permettant d'afficher les instantannes.
% Necessite en entree:
%         1. matrices coordonnees: X et Y 
%         2. valeurs: 
%         3. les valeurs extremes a considerer (pour avoir un affichage a la
%         meme echelle de couleur sur les differents instantannes
%         4. on definit les axes: [xmin xmax]x[ymin ymax]

X = [X;X];
Y = [-Y(end:-1:1,:);Y];
M = [M(end:-1:1,:);M];
% keyboard;
if (nargin==3)
  cmin=min(min(M));
  cmax=max(max(M));
  xmin=min(min(X));
  ymin=min(min(Y));
  xmax=max(max(X));
  ymax=max(max(Y));
elseif (nargin==5)
  xmin=min(min(X));
  ymin=min(min(Y));
  xmax=max(max(X));
  ymax=max(max(Y));
end
format long
%hold on;
h=figure;
set(h,'paperposition', [0.25 2.5 5.675 4.95])
%h = subplot('position',[.05 .55 .85 .4]);
%h = subplot('position',[.05 .05 .85 .4]);
hh=pcolor(X,Y,M);
shading interp;
caxis([cmin cmax])
colormap(jet);
colorbar
axis ([xmin xmax ymin ymax])
axis image
drawnow
