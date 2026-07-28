function film2D(base, ext, debut, fin, cmin, cmax, dessine_contour, num_inst, raff)

if (nargin <= 8)
  raff = 1;
end
if (nargin <= 7)
  num_inst = 1;
end
if (nargin <= 6)
  dessine_contour = 0;
end

h=gcf;%figure(1);
set(h,'visible','on');
set(h,'paperposition', [0.25 2.5 5.675 4.95]);
set(h,'paperpositionmode','auto');
%set(1,'Position',[80 80 800 400]);
mode = '-djpeg90';

% on charge le dernier
if (dessine_contour == 1)
  txt = [base, entier_to_string(fin), ext, '.dat'];
  [XI,YI,Z,coor,V] = loadND(txt, num_inst, raff);
  [non_loc, num_inter] = contour2D(XI,YI,V);
end
%keyboard;

for i = debut:fin
  entier = entier_to_string(i);
  fic = strcat(base, entier); fic = strcat(fic,ext);
  fic = strcat(fic,'.dat')
  [XI,YI,ZI,coor,V] = loadND(fic, num_inst, raff);
  if (dessine_contour == 1)
    V(non_loc) = NaN;
    V(num_inter) = NaN;
  end
  fics = strcat(base, entier); fics = strcat(fics,ext);
  fics = strcat(fics,'.jpg')
  plot2dinst(XI, YI, V, cmin, cmax, 1);
  print(h, mode, fics);
  %close(1)
  % pause
end

