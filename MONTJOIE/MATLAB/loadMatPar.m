function A = loadMat(fic, n, assemble)

if (nargin == 2)
  assemble = 1;
end

m = 0;
for i=1:n
  A = loadMat([fic, '_P', num2str(i-1), '.dat']);
  m = max(m, size(A, 1));
  m = max(m, size(A, 2));
end

A = sparse(m, m);
for i=1:n
  B = loadMat([fic, '_P', num2str(i-1), '.dat']);
  if (assemble == 1)
     A(1:size(B,1), 1:size(B,2)) = A(1:size(B,1), 1:size(B,2)) + B;
  else
     A(1:size(B,1), 1:size(B,2)) = B;
  end
end
