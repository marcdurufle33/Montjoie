function X=load1D(fic)
  X = [];
  fid = fopen(fic,'r');
  if (fid ~= -1)
    n = fread(fid,1,'int');
    data = fread(fid,2*n,'double');
    X = data(1:2:(2*n))+1i*data(2:2:(2*n));
    fclose(fid);
  end