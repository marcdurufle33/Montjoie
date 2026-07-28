function [sor] = plot3d_plane(X,Y,Z,M,coor,cmin,cmax,num_figure)

sor=[];

% Fonction permettant d'afficher les instantannes.
% Necessite en entree:
%         1. matrices coordonnees: X, Y et Z
%         2. valeurs: M
%         3. coordonnee du plan de coupe
%         4. les valeurs extremes a considerer (pour avoir un affichage a la
%         meme echelle de couleur sur les differents instantannes
%         5. on definit les axes: [xmin xmax]x[ymin ymax]x[zmin
%         zmax]

if (nargin==5)
  cmin=min(min(min(M)));
  cmax=max(max(max(M)));
  xmin=min(min(min(X)));
  xmax=max(max(max(X)));
  ymin=min(min(min(Y)));
  ymax=max(max(max(Y)));
  zmin=min(min(min(Z)));
  zmax=max(max(max(Z)));
  num_figure = 0;
elseif (nargin>=7)
  xmin=min(min(min(X)));
  xmax=max(max(max(X)));
  ymin=min(min(min(Y)));
  ymax=max(max(max(Y)));
  zmin=min(min(min(Z)));
  zmax=max(max(max(Z)));
  if (nargin <= 7)
    num_figure = 0;
  end
end

if (num_figure == 0)
  h = figure;
  set(h,'paperposition', [0.25 2.5 5.675 4.95]);
  set(h,'color',[1 1 1]);
else
  clf(num_figure);
end

subplot(2,2,1);
hold on;

subplot(2,2,2);
Val = M(:,:,1);
%keyboard;
pcolor(linspace(ymin,ymax,size(Val,1)),linspace(zmin,zmax,size(Val,2)),transpose(Val)); 
caxis([cmin cmax]); colormap(jet(20000));
axis ([ymin ymax zmin zmax]); shading interp;
axis image;

subplot(2,2,3);
Val = M(:,:,2);
pcolor(linspace(xmin,xmax,size(Val,2)),linspace(zmin,zmax,size(Val,1)),Val);
caxis([cmin cmax]); colormap(jet(20000));
axis ([xmin xmax zmin zmax]); shading interp;
axis image;

subplot(2,2,4);
Val = M(:,:,3);
pcolor(linspace(xmin,xmax,size(Val,2)),linspace(ymin,ymax,size(Val,1)),Val);
caxis([cmin cmax]); colormap(jet(20000));
axis ([xmin xmax ymin ymax]); shading interp;
axis image;