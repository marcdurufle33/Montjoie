function A = load_full(fic)
  fid = fopen(fic,'r');
  m = fread(fid,1,'int');
  n = fread(fid,1,'int');
  data = fread(fid,2*m*n,'double');
  X = data(1:2:(2*m*n)) + 1i*data(2:2:(2*m*n));
  X = transpose(X);
  fclose(fid);
  A = reshape(X,n,m); A = transpose(A);
  