function [non_loc, num_inter] = contour2D(XI,YI,V)

non_loc = []; num_inter = []; seuil = 1e-80;
% on elime tous les ilots
nx = length(XI);
for i = 2: (nx-1)
  for j = 2: (nx-1)
    if ((abs(V(i+1,j)) < seuil )&&(abs(V(i-1,j)) < seuil)&&(abs(V(i,j+1)) < seuil)&&(abs(V(i,j-1)) < seuil))
        V(i,j) = 0;
    end
  end
end
% on cherche les points qui n'ont pas ete localises
non_loc = find(abs(V)<seuil);
% tous les points a l'interface sont marques
nb = 0;
for i = 2: (nx-1)
  for j = 2: (nx-1)
    if (abs(V(i,j)) < seuil)
      k = (j-1)*nx + i;
      interface = 0;
      if (abs(V(i+1,j)) > seuil)
        interface = 1;
      end
      if (abs(V(i-1,j)) > seuil)
        interface = 1;
      end
      if (abs(V(i,j+1)) > seuil)
        interface = 1;
      end
      if (abs(V(i,j-1)) > seuil)
        interface = 1;
      end
      if (interface == 1)
        nb = nb + 1;
        num_inter(nb) = k;
      end
    end
  end
end
%keyboard;
