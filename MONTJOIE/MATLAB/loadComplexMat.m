function A = loadComplexMat(fic)
fid = fopen(fic);
A = fscanf(fid,'%d %d (%g,%g)');
B = transpose(reshape(A, 4, length(A)/4));
C(:,1:2) = B(:,1:2);
C(:, 3) = B(:,3) + 1i*B(:,4);
A = spconvert(C);
