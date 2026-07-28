function [X,Y,Z,coor,Vx,Vy,Vz,dVx,dVy,dVz] = loadE(fic, num_inst, raff)

  % if the grid number is not specified, we take the first one
  % if the refinement ratio is not specified, we take one
  if (nargin == 1)
    num_inst = 1;
    raff = 1;
  end
  if (nargin == 2)
    raff = 1;
  end
  
  [X,Y,Z,coor,Vx] = loadND([fic,'_U0.dat'], num_inst, raff);
  [X,Y,Z,coor,Vy] = loadND([fic,'_U1.dat'], num_inst, raff);
  [X,Y,Z,coor,Vz] = loadND([fic,'_U2.dat'], num_inst, raff);
  [X,Y,Z,coor,dVx] = loadND([fic,'_dU0.dat'], num_inst, raff);
  [X,Y,Z,coor,dVy] = loadND([fic,'_dU1.dat'], num_inst, raff);
  [X,Y,Z,coor,dVz] = loadND([fic,'_dU2.dat'], num_inst, raff);
  