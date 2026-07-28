function vitesse=calcule_vitesse(vect,Deltat)

vect = vect(:)';
vect2=circshift(vect,[0 2]);
%vect2=circshift(vect,[0 1]);
vitesse=(vect-vect2);
vitesse=vitesse(3:length(vitesse));
vitesse=vitesse/Deltat/2;
