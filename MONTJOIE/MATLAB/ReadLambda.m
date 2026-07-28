function [L, tps, freq, amo] = ReadLambda(name_file,Nmax,MAXI)

%[s,m]=system(['ls -l ' name_file ]);
%for i=1:5
%    [b,m]=strtok(m,' ');
%end

fid = fopen(name_file);

nev = fread(fid, 1, 'int');
freq = fread(fid, nev, 'float');

nev = fread(fid, 1, 'int');
amo = fread(fid, nev, 'float');

if (Nmax > nev)
  Nmax = nev)
end

L = zeros(MAXI, Nmax);
tps = zeros(MAXI, 1)

compt=0;
try
while(~feof(fid) && compt<MAXI)
   compt = compt+1;
   n = fread(fid, 1, 'int');
   Lambda = fread(fid, n, 'float');
   tps(compt) = Lambda(1);
   L(compt,:) = Lambda(2:Nmax+1);
   
end
end

L = L(1:compt-1, 1:Nmax);
tps = tps(1:compt-1);
