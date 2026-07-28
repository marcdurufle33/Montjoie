function [Points, L] = ReadFarField(name_file, MAXI)

fid = fopen(name_file);

% type of data (0 : float, 1 : double, 2 : complex float, 3 : complex double)
type_data = fread(fid, 1, 'int');
prec = 'double'; prec_size = 8;
if ((type_data == 0)||(type_data == 2))
  prec = 'float'; prec_size = 4;
end

n = fread(fid, 1, 'int'); n=min(n,3e7);
Points = fread(fid, n, prec);
Points = reshape(Points, 3, n/3);

L = zeros(MAXI, n/3+1);

compt=0;
try
  while(~feof(fid) && compt<MAXI)
     compt=compt+1;
     n = fread(fid, 1, 'int');
     Vloc = fread(fid, n, prec);
     L(compt, :) = Vloc;   
  end
end

L = L(1:compt-1, :);

V = L;
[a, b] = mon_unique(V(:,1));
L = zeros(length(a), size(V, 2));
L(:,1) = a;
for i = 2:size(V,2)
  L(:,i) = V(b,i);
end

