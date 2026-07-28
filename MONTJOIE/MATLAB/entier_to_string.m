function entier = entier_to_string(i)


fchaine4 = num2str(mod(i,10));
n0 = floor(i/10.0);
fchaine3 = num2str(mod(n0,10));
n0 = floor(n0/10.0);
fchaine2 = num2str(mod(n0,10));
n0 = floor(n0/10.0);
fchaine1 = num2str(mod(n0,10));
n0 = floor(n0/10.0);
entier = strcat(fchaine1,fchaine2,fchaine3,fchaine4);

