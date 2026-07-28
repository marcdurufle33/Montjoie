function sor = write_full_s(fic, nu)
  fid = fopen(fic,'w');
  s = size(nu);
  % writing matrix in binary format (single precision)
  fwrite(fid, s(1), 'int');
  fwrite(fid, s(2), 'int');
  M = reshape(transpose(nu), s(1)*s(2), 1);
  fwrite(fid, M, 'float');
  fclose(fid);
  
  