function [sor] = plot3d_piano_defil(base,instant,coef,zplaque,cmin,cmax,num_figure)
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
file1= [base '/SteinwayCoupe' entier_to_string(instant) '_U0.dat']
file2 = [base '/SteinwayPlate' entier_to_string(instant) '_U2.dat']
file3 = [base '/../Sismo_Last_Point0.txt'];
ecriture = [base '/Slice3D' entier_to_string(instant)];
if(~exist(file1))
    'je t insulte'
    return
end
if(~exist(file2))
    'je t insulte deux'
    return
end
try
   sismo = load(file3);
catch
   'je t insulte trois'
    return
end



%% en bas

[X,Y,Z,coor,M] = loadND(file1);
[Xp,Yp,Zp,coorp,Vp]=loadND(file2);
  [non_loc, num_inter] = contour2D(Xp,Yp,Vp);
    Vp(non_loc) = NaN;
  %  Vp(num_inter) =1e9;

zplaque_a_definir=0;
if (nargin==2)
  cmin=min(min(min(M)));
  cmax=max(max(max(M)));
  xmin=min(min(min(X)));
  xmax=max(max(max(X)));
  ymin=min(min(min(Y)));
  ymax=max(max(max(Y)));
  zmin=min(min(min(Z)));
  zmax=max(max(max(Z)));
  num_figure = 0;
  coef=1e3;
  zplaque_a_definir=1;
elseif (nargin==3)
  cmin=min(min(min(M)));
  cmax=max(max(max(M)));
  xmin=min(min(min(X)));
  xmax=max(max(max(X)));
  ymin=min(min(min(Y)));
  ymax=max(max(max(Y)));
  zmin=min(min(min(Z)));
  zmax=max(max(max(Z)));
  num_figure = 0;
  zplaque_a_definir=1;
elseif (nargin==4)
  cmin=min(min(min(M)));
  cmax=max(max(max(M)));
  xmin=min(min(min(X)));
  xmax=max(max(max(X)));
  ymin=min(min(min(Y)));
  ymax=max(max(max(Y)));
  zmin=min(min(min(Z)));
  zmax=max(max(max(Z)));
  num_figure = 0;

elseif (nargin>=6)
  xmin=min(min(min(X)));
  xmax=max(max(max(X)));
  ymin=min(min(min(Y)));
  ymax=max(max(max(Y)));
  zmin=min(min(min(Z)));
  zmax=max(max(max(Z)));
  if (nargin == 6)
    num_figure = 0;
  
  end
end
cmin = -max(abs(cmin),abs(cmax));
cmax = -cmin;
sor=[sor,cmin,cmax];
nb_div = 2;
%keyboard;
[val,i0] = min(abs(X-coor(1)));
[val,j0] = min(abs(Y-coor(2)));
[val,k0] = min(abs(Z-coor(3)));
i0b = length(i0:-nb_div:1);
j0b = length(j0:-nb_div:1);
k0b = length(k0:-nb_div:1);
ix = [fliplr(i0:-nb_div:1),(i0+nb_div:nb_div:length(X))];
iy = [fliplr(j0:-nb_div:1),(j0+nb_div:nb_div:length(Y))];
iz = [fliplr(k0:-nb_div:1),(k0+nb_div:nb_div:length(Z))];
xc = X(ix);
yc = Y(iy);
zc = Z(iz);
%[x,y,z] = ndgrid(xc, yc, zc);
W = zeros(length(xc),length(yc),length(zc));
W(:,i0b,:) = M(iy,iz,1);
W(j0b,:,:) = transpose(M(ix,iz,2));
W(:,:,k0b) = M(ix,iy,3);
if (num_figure == 0)
  h = figure(1); clf;
  %set(h,'paperposition', [0.25 2.5 5.675 4.95]);
  set(h,'color',[1 1 1]);
  
else
  clf(num_figure);
end




%% en bas
subplot(2,1,2);
%keyboard;
if(zplaque_a_definir)
    zplaque = Z(k0);
end

slice_piano(xc,yc,zc,W, X(i0),Y(j0),Z(k0),Xp,Yp,Zp,coef*Vp,zplaque);%,'nearest');
shading interp;
hold on;
axis image;


% traits noirs
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
plot3(X1,Y1,[zplaque,zplaque],'k')

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
plot3(X1,Y1,[zplaque,zplaque],'k')


caxis([cmin cmax]);
colormap(jet);
% colorbar;
%axis ([xmin xmax ymin ymax zmin zmax]);
axis ([xmin xmax ymin ymax zmin zmax]);
zoom(1.75);
axis off;
grid off;
%view(300,45);
view(320,35);
%view(200,-20);
drawnow;

%snam='PowerPoint';
%s=hgexport('readstyle',snam);
%hgexport(gcf,[ file1 '.eps' ],s);


%% en haut
% pas de temps entre les interpolations
Deltat = 1e-5;

tbegin = 0;
tend = min(3e-2,sismo(end,1));

[mi,indbegin]=min((sismo(:,1)-tbegin).^2);
[mi,indend]=min((sismo(:,1)-tend).^2);
indend = indend-4;
Deltatsismo = 1e-6;%sismo(2,1);
acce =  calcule_vitesse(calcule_vitesse(sismo(:,2)',Deltatsismo),Deltatsismo);

h1=subplot(2,1,1); cla(h1); axes(h1);
plot(sismo(indbegin:indend,1),acce(indbegin:indend),'-b.',...
                'LineWidth',2,...
                'MarkerEdgeColor','b',...
                'MarkerSize',0.1); %sismo(indbegin:indend,2)
xlim([tbegin*0.99 tend*1.01]);
ylabel('Bridge acceleration (m/s^2)');

M=max(abs(acce(indbegin:indend)));%sismo(indbegin:indend,2)));%
hold on, plot([instant*Deltat,instant*Deltat],[-M,M],'Color','r','LineWidth',3);
title(['t = ' num2str(instant*Deltat,'%2.5f')]);


%%  impression
%mode1 = '-djpeg90'
%print(gcf, mode1, ecriture);
% exportpspdf(gcf,ecriture);
