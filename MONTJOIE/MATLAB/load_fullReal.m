function A = load_fullReal(fic)
  fid = fopen(fic,'r');
  m = fread(fid,1,'int');
  n = fread(fid,1,'int');
  data = fread(fid,m*n,'double');
  X = transpose(data);
  fclose(fid);
  A = reshape(X,n,m); A = transpose(A);
  