function X = load1D_float(fic)
  fid = fopen(fic,'r');
  n = fread(fid,1,'int');
  X = fread(fid,n,'float');
  fclose(fid);
