function A = affiche_mode(root, ext, n1, n2)
  place=strfind(root,'/');
  DOSSIER = root(1:place(end));
for i = n1:n2
  [root, entier_to_string(i), ext, '.dat']
  [X,Y,Z,coor,V] = loadND([root, entier_to_string(i), ext, '.dat']);
  if (numel(Z) < 2)
    plot2dinst(X, Y, real(V));
  else
    plot3d_planec(X, Y, Z, real(V), coor);
  end
  pause
end
