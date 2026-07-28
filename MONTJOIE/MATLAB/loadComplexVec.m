function A = loadComplexVec(fic)
fid = fopen(fic);
A = fscanf(fid,'(%g,%g)');

