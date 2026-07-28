function [non_loc, num_inter] = find_contour(XI,YI,ZI,coor,V)

non_loc = []; num_inter = []; seuil = 1e-20;
% on elime tous les ilots
nx = length(XI);
seuil = 1e-20;
for i = 2: (nx-1)
  for j = 2: (nx-1)
    for num = 1:3
      if ((abs(V(i+1,j,num)) < seuil )&&(abs(V(i-1,j,num)) < seuil)&&(abs(V(i,j+1,num)) < seuil)&&(abs(V(i,j-1,num)) < seuil))
        V(i,j,num) = 0;
      end
    end
  end
end
% on cherche les points qui n'ont pas ete localises
non_loc = find(abs(V)<seuil);
% tous les points a l'interface sont marques
nb = 0;
for i = 2: (nx-1)
  for j = 2: (nx-1)
    for num = 1:3
      if (abs(V(i,j,num)) < seuil)
        if ((num == 2)&&(j<=30)&&(i<=180)&&(i>=150))
        else
        k = (num-1)*nx*nx + (j-1)*nx + i;
        interface = 0;
        if (abs(V(i+1,j,num)) > seuil)
          interface = 1;
        end
        if (abs(V(i-1,j,num)) > seuil)
          interface = 1;
        end
        if (abs(V(i,j+1,num)) > seuil)
          interface = 1;
        end
        if (abs(V(i,j-1,num)) > seuil)
          interface = 1;
        end
        if (interface == 1)
          nb = nb + 1;
          num_inter(nb) = k;
        end
      end
      end
    end
  end
end
%keyboard;