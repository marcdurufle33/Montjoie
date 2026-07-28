function U = loadSismo(nom_fichier)

V = load(nom_fichier);
[a, b] = mon_unique(V(:,1));
U = zeros(length(a), size(V, 2));
U(:,1) = a;
for i = 2:size(V,2)
  U(:,i) = V(b,i);
end

