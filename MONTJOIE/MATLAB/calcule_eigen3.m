% on charge les matrices
L0 = loadComplexMat('../mat_jacobi_P0.dat');
L1 = loadComplexMat('../mat_jacobi_P1.dat');
L2 = loadComplexMat('../mat_jacobi_P2.dat');

% les numeros pour Omega 1
dof0 = load('../num_gamma0_P1.dat')+1;
dof1 = load('../num_gamma1_P1.dat')+1;

% on recupere les blocs
A10 = full(L0);
A01 = full(L1(dof0, dof0));
A02 = full(L1(dof0, dof1));
A31 = full(L1(dof1, dof0));
A32 = full(L1(dof1, dof1));
A23 = full(L2);

clear L0;
clear L1;
clear L2;

M = [A01*A10, A02*A23; A31*A10, A32*A23];
lambda = sqrt(eig(M));
 