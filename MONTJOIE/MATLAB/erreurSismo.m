function err = erreurSismo(U, V)
 
for nb = 2:size(V,2)
  Ui = interp1(U(:,1),U(:,nb),V(:,1),'cubic');
  err(nb) = sum((Ui-V(:,nb)).^2);
  absolu(nb) = sum((V(:,nb)).^2);
end

err = sqrt(sum(err))/sqrt(sum(absolu));
