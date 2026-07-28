function film3D(base, ext, debut, fin, cmin, cmax, dessine_contour, num_inst, raff)

if (nargin == 8)
  raff = 1;
end
if (nargin == 7)
  raff = 1;
  num_inst = 1;
end
if (nargin == 6)
  dessine_contour = 0;
  raff = 1;
  num_inst = 1;
end
figure(1);
set(1,'visible','off');
set(1,'paperposition', [0.25 2.5 5.675 4.95]);
set(1,'paperpositionmode','auto');
%set(1,'Position',[80 80 800 400]);
mode = '-djpeg90';

% on charge le dernier
if (dessine_contour == 1)
  txt = [base, entier_to_string(fin), ext, '.dat'];
  [XI,YI,ZI,coor,V] = loadND(txt);
  [non_loc, num_inter] = find_contour(XI,YI,ZI,coor,V);
end

for i = debut:fin
  entier = entier_to_string(i);
  fic = strcat(base, entier); fic = strcat(fic,ext);
  fic = strcat(fic,'.dat')
  [XI,YI,ZI,coor,V] = loadND(fic, num_inst, raff);
  if (dessine_contour == 1)
    V(non_loc) = 0;
    V(num_inter) = NaN;
  end
  fics = strcat(base, entier);  fics = strcat(fics,ext);
  fics = strcat(fics,'.jpg')
  plot3d_planec(XI, YI, ZI, real(V), coor, cmin, cmax, 1);
  print(1, mode, fics);
  %close(1)
end

