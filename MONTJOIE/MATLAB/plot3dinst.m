function sor=plot3dinst(X,Y,Z,M,coor,cmin,cmax,xmin,xmax,ymin,ymax,zmin,zmax);
% Fonction permettant d'afficher les instantannes.
% Necessite en entree:
%         1. matrices coordonnees: X, Y et Z
%         2. valeurs: M
%         3. coordonnee du plan de coupe
%         4. les valeurs extremes a considerer (pour avoir un affichage a la
%         meme echelle de couleur sur les differents instantannes
%         5. on definit les axes: [xmin xmax]x[ymin ymax]x[zmin
%         zmax]
format long
if (nargin==5)
  %keyboard;
  cmin=min(min(min(M)));
  cmax=max(max(max(M)));
  xmin=min(min(min(X)));
  xmax=max(max(max(X)));
  ymin=min(min(min(Y)));
  ymax=max(max(max(Y)));
  zmin=min(min(min(Z)));
  zmax=max(max(max(Z)));
elseif (nargin==7)
  xmin=min(min(min(X)));
  xmax=max(max(max(X)));
  ymin=min(min(min(Y)));
  ymax=max(max(max(Y)));
  zmin=min(min(min(Z)));
  zmax=max(max(max(Z)));
end

h = figure;
set(h,'paperposition', [0.25 2.5 5.675 4.95]);
set(h,'color',[1 1 1])
subplot(2,2,1)
hh=slice(Y,X,Z,M,[coor(1) coor(1)],coor(2),[coor(3) coor(3)]);
shading interp;
hold on;
axis image;
X1(1)=coor(1);X1(2)=coor(1);
Y1(1)=coor(2);Y1(2)=coor(2);
Z1(1)=zmin;Z1(2)=zmax;
plot3(X1,Y1,Z1,'k');
hold on;
X1(1)=xmin;X1(2)=xmin;
Y1(1)=coor(2);Y1(2)=coor(2);
Z1(1)=zmin;Z1(2)=zmax;
plot3(X1,Y1,Z1,'k');
hold on;
X1(1)=xmax;X1(2)=xmax;
Y1(1)=coor(2);Y1(2)=coor(2);
Z1(1)=zmin;Z1(2)=zmax;
plot3(X1,Y1,Z1,'k');
hold on;
X1(1)=coor(1);X1(2)=coor(1);
Y1(1)=ymin;Y1(2)=ymin;
Z1(1)=zmin;Z1(2)=zmax;
plot3(X1,Y1,Z1,'k');
hold on;
X1(1)=coor(1);X1(2)=coor(1);
Y1(1)=ymax;Y1(2)=ymax;
Z1(1)=zmin;Z1(2)=zmax;
plot3(X1,Y1,Z1,'k');
hold on;
X1(1)=coor(1);X1(2)=coor(1);
Y1(1)=ymin;Y1(2)=ymax;
Z1(1)=zmin;Z1(2)=zmin;
plot3(X1,Y1,Z1,'k');
hold on;
X1(1)=coor(1);X1(2)=coor(1);
Y1(1)=ymin;Y1(2)=ymax;
Z1(1)=zmax;Z1(2)=zmax;
plot3(X1,Y1,Z1,'k');
hold on;
X1(1)=xmin;X1(2)=xmax;
Y1(1)=coor(2);Y1(2)=coor(2);
Z1(1)=zmax;Z1(2)=zmax;
plot3(X1,Y1,Z1,'k');
hold on;
X1(1)=xmin;X1(2)=xmax;
Y1(1)=coor(2);Y1(2)=coor(2);
Z1(1)=zmin;Z1(2)=zmin;
plot3(X1,Y1,Z1,'w');
hold on;
caxis([cmin cmax]);
colormap(jet);
% colorbar;
axis ([xmin xmax ymin ymax zmin zmax]);
axis off;
grid off;
view(290,45)
drawnow
%keyboard
subplot(2,2,3)
hh=slice(X,Y,Z,M,[coor(1) coor(1) coor(1)],coor(2),[coor(3) coor(3)]);
shading interp;
hold on;
%axis equal;
axis image;
caxis([cmin cmax]);
colormap(jet);
% colorbar;
axis ([xmin xmax ymin ymax zmin zmax]);
axis off;
grid off;
view(270,0)
drawnow
subplot(2,2,2)
hh=slice(X,Y,Z,M,[coor(1) coor(1) coor(1)],coor(2),[coor(3) coor(3)]);
shading interp;
hold on;
axis equal;
axis image;
% plot3(x,y,z,'w','linewidth',2)
% plot3(x/2,y/2,z,'w','linewidth',2)
caxis([cmin cmax]);
colormap(jet);
% colorbar;
axis ([xmin xmax ymin ymax zmin zmax]);
axis off;
grid off;
view(270,90)
drawnow
subplot(2,2,4)
hh=slice(X,Y,Z,M,[coor(1) coor(1) coor(1)],coor(2),[coor(3) coor(3)]);
shading interp;
hold on;
axis equal;
axis image;
caxis([cmin cmax]);
colormap(jet);
% colorbar;
axis ([xmin xmax ymin ymax zmin zmax]);
axis off;
grid off;
view(0,0)
drawnow

