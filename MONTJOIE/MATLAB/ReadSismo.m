function U = ReadSismo(name_file, MAXI, no_time)

if (nargin == 1)
  MAXI = 10000000;
  no_time = false;
end

if (nargin == 2)
  no_time = false;
end

fid = fopen(name_file);

% type of data (0 : float, 1 : double)
type_data = fread(fid, 1, 'int');
prec = 'double';
if (type_data == 0)
  prec = 'float';
end

compt=0;
try
while(~feof(fid) && compt<MAXI)
   compt = compt+1;
   n = fread(fid, 1, 'int');
   V = fread(fid, n, prec);
   if (compt == 1)
     sis = zeros(MAXI, n);
   end
   sis(compt,:) = V;
end
end

if (no_time == false)
  sis = sis(1:compt-1, :);
  [a, b] = mon_unique(sis(:,1));
  U = zeros(length(a), size(sis, 2));
  U(:,1) = a;
  for i = 2:size(sis,2)
    U(:,i) = sis(b,i);
  end
else
  U = sis(1:compt-1, :);
end

