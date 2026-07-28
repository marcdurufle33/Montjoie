function sor = write_full_c(fic, nu)
  fid = fopen(fic,'w');
  s = size(nu);
  % writing matrix in binary format (double precision)
  fwrite(fid, s(1), 'int');
  fwrite(fid, s(2), 'int');
  M = reshape(transpose(nu), s(1)*s(2),1);
  Mr = zeros(s(1)*s(2)*2, 1);
  Mr(1:2:end) = real(M);
  Mr(2:2:end) = imag(M);
  fwrite(fid, Mr, 'float');
  fclose(fid);
  
  
