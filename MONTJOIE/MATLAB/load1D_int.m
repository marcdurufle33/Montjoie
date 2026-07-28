function X=load1D(fic)
  fid=fopen(fic,'r');
  n=fread(fid,1,'int');
  X = fread(fid,n,'int');
  fclose(fid);
