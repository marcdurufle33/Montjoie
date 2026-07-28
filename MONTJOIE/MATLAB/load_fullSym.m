function A = load_fullSym(fic)
  fid = fopen(fic,'r');
  n = fread(fid,1,'int');
  n = fread(fid,1,'int');
  data = fread(fid,n*(n+1),'double');
  X = data(1:2:((n+1)*n)) + 1i*data(2:2:((n+1)*n));
  X = transpose(X);
  fclose(fid);
  A = zeros(n,n);
  ind = 1;
  for i=1:n
    taille = n-i;
    A(i,i:n) = X(ind:(ind+taille) ).';
    A(i:n,i) = X(ind:(ind+taille) );
    ind = ind + taille + 1;
  end
  