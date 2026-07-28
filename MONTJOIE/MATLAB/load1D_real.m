function X = load1D_real(fic)
  fid=fopen(fic,'r');
  n=fread(fid,1,'int');
  X = fread(fid,n,'double');
  fclose(fid);
