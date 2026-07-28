#!/usr/bin/env python 
""" Module for reading Montjoie output files, and displaying them,
   Python equivalent of Matlab scripts"""

import numpy as np
import scipy
import scipy.sparse
import os, sys, pickle
import pylab
import struct, wave
from pylab import cos, sin, pi
from scipy.interpolate import RectBivariateSpline
from scipy.interpolate import interp1d
import matplotlib.colors

def Refine1D(xmin, xmax, U, raff):
  """ interpolates a 1-D field with cubic interpolation 
  usage : V = Refine1D(xmin, xmax, U, raff)
         raff : level of refinement """
  nbx = len(U)
  x = pylab.linspace(xmin, xmax, nbx);
  xr = pylab.linspace(xmin, xmax, raff*(nbx-1)+1);
  Up = interp1d(x, U, 'cubic')
  return Up(xr)

def Refine2D(xmin, xmax, ymin, ymax, U, raff):
  """ interpolates a 2-D field with cubic interpolation 
    usage : V = Refine2D(xmin, xmax, ymin, ymax, U, raff)
         raff : level of refinement """
  nbx = U.shape[0]
  nby = U.shape[1]
  x = pylab.linspace(xmin, xmax, nbx);
  y = pylab.linspace(ymin, ymax, nby);
  xr = pylab.linspace(xmin, xmax, raff*(nbx-1)+1);
  yr = pylab.linspace(ymin, ymax, raff*(nby-1)+1);
  if (U.dtype == 'complex128'):          
    ip_r = RectBivariateSpline(x, y, U.real)
    ip_i = RectBivariateSpline(x, y, U.imag)
    U = ip_r(xr, yr) + 1j*ip_i(xr, yr)
  else:
    ip = RectBivariateSpline(x, y, U)
    U = ip(xr, yr)

  return U

def RefineFft(u, s):
  """ Interpolates a field u to fit a new vector of size s, it works for smooth period field only
  usage : yinterp = RefineFft(y, N)
  
  y : periodic data of size N such as y[0] = y[N-1]
      N-1 should be taken as a factor of powers of 2, 3, 5 or 7 for best efficiency
  s : new size (again s-1 should be taken as a factor of powers of 2, etc) """
  
  
  N = len(u)
  # we perform fft by excluding the last point (by periodicity)
  uchap = pylab.fft(u[0:N-1])
  uchap_pad = pylab.zeros(s-1) + 1j*pylab.zeros(s-1)
  # we fill uchap_pad with additional zeros (for high frequency)
  uchap_pad[0:(N-1)//2] = uchap[0:(N-1)//2]
  uchap_pad[s-1-(N-1)//2:] = uchap[(N-1)//2:]
  ur = pylab.zeros(s) + 1j*pylab.zeros(s)
  # we perform inverse fft to obtain ur on refined points
  ur[0:s-1] = (s-1) / float(N-1) * pylab.ifft(uchap_pad)
  # periodic condition
  ur[s-1] = ur[0]
  # we return only real part (u is assumed to be a real vector)
  return pylab.real(ur)

def loadND(nom_fichier, num_inst=1, raff=1):
  """ Reading of a "Montjoie" output file (MATLAB output) 
  usage : [X, Y, Z, coor, V] = loadND(nom_fichier)
          [X, Y, Z, coor, V] = loadND(nom_fichier, num_inst)
          [X, Y, Z, coor, V] = loadND(nom_fichier, num_inst, raff)

  nom_fichier : name of the output file
  num_inst : number of the snapshot (1 if not specified)
  raff : level of refinement
  X : x-coordinates
  Y : y-coordinates
  Z : z-coordinates
  coor : center of three planes (for FileOutputGrille)
  V : values of the solution on points (x,y,z) """ 
  
  # initialisation of output arrays
  XI = []; YI = []; ZI = []; coor = []; V = [];
  
  # we read the number of grids 
  fileobj = open(nom_fichier, mode='rb')
  nb_grids = pylab.fromfile(fileobj, 'i', 1)[0]
  
  # the dimension
  dim = pylab.fromfile(fileobj, 'i', 1)[0]
  
  # type of data (0 : float, 1 : double, 2 : complex float, 3 : complex double)
  type_data = pylab.fromfile(fileobj, 'i', 1)[0]
  prec = 'd'; prec_size = 8
  if ((type_data == 0) or (type_data == 2)):
    prec = 'f';
    prec_size = 4
  
  # cplx = 2 for complex and 1 for real numbers
  cplx = 1;
  if (type_data >= 2):
    cplx = 2;
  
  prec_size = prec_size*cplx;
  
  # type of grid
  type_grid = pylab.fromfile(fileobj, 'i', 1)
  
  if (dim == 1):
    xmin = pylab.fromfile(fileobj, prec, 1)[0]
    xmax = pylab.fromfile(fileobj, prec, 1)[0]
    nbx = pylab.fromfile(fileobj, 'i', 1)[0]
    n = pylab.fromfile(fileobj, 'i', 1)[0]
    taille = int(cplx*nbx);
    data = pylab.fromfile(fileobj, prec, taille)
    if (cplx == 2):
      M = data[0:taille:2] + 1j*data[1:taille:2];
    else:
      M = data;
      
    fileobj.close()
    V = M
    if (raff != 1):
      V = Refine1D(xmin, xmax, V, raff)
      
    XI = pylab.linspace(xmin, xmax, (nbx-1)*raff + 1)
  elif (dim == 2):
    
    # depending the type of grid, we read the appropriate datas
    if (type_grid == 0):

      # SismoPlane
      offset = 0;
      xmin_ = pylab.zeros(nb_grids); xmax_ = pylab.zeros(nb_grids);
      ymin_ = pylab.zeros(nb_grids); ymax_ = pylab.zeros(nb_grids);
      nbx_ = pylab.zeros(nb_grids).astype('int'); nby_ = pylab.zeros(nb_grids).astype('int');
      for i in range(nb_grids):
        xmin_[i] = pylab.fromfile(fileobj, prec, 1)
        xmax_[i] = pylab.fromfile(fileobj, prec, 1)
        nbx_[i] = pylab.fromfile(fileobj, 'i', 1)[0]
        ymin_[i] = pylab.fromfile(fileobj, prec, 1)
        ymax_[i] = pylab.fromfile(fileobj, prec, 1)
        nby_[i] = pylab.fromfile(fileobj, 'i', 1)[0]     
        if (i < num_inst-1):
          offset = offset + nbx_[i]*nby_[i];
        
       
      xmin = xmin_[num_inst-1]; xmax = xmax_[num_inst-1];
      ymin = ymin_[num_inst-1]; ymax = ymax_[num_inst-1];
      nbx = nbx_[num_inst-1]; nby = nby_[num_inst-1];
      n = pylab.fromfile(fileobj, 'i', 1)[0]
      taille = int(cplx*nbx*nby);
      # we read values related to the required snapshot
      fileobj.seek(offset*prec_size, 1);
      data = pylab.fromfile(fileobj, prec, taille)
      if (cplx == 2):
        M = data[0:taille:2] + 1j*data[1:taille:2];
      else:
        M = data;
      
      fileobj.close()

      # here we obtain the transpose of the data written by the user
      # it is because 2-D displays naturally transpose the matrix (rows are displayed along y)
      U = pylab.reshape(M, (nby, nbx));
      if (raff != 1):
        U = Refine2D(ymin, ymax, xmin, xmax, U, raff)
      
      X = pylab.linspace(xmin, xmax, (nbx-1)*raff+1);
      Y = pylab.linspace(ymin, ymax, (nby-1)*raff+1);
      V = U;
      XI = X; YI = Y;
              
    elif (type_grid == 1):
      
      # SismoLine
      offset = 0;
      xmin_ = pylab.zeros(nb_grids); xmax_ = pylab.zeros(nb_grids);
      ymin_ = pylab.zeros(nb_grids); ymax_ = pylab.zeros(nb_grids);
      nbx_ = pylab.zeros(nb_grids).astype('int');
      for i in range(nb_grids):
        xmin_[i] = pylab.fromfile(fileobj, prec, 1)[0]
        xmax_[i] = pylab.fromfile(fileobj, prec, 1)[0]
        ymin_[i] = pylab.fromfile(fileobj, prec, 1)[0]
        ymax_[i] = pylab.fromfile(fileobj, prec, 1)[0]
        nbx_[i] = pylab.fromfile(fileobj, 'i', 1)[0]
        if (i < num_inst-1):
          offset = offset + nbx_[i];
        
      
      xmin = xmin_[num_inst-1]; xmax = xmax_[num_inst-1];
      ymin = ymin_[num_inst-1]; ymax = ymax_[num_inst-1];
      nbx = nbx_[num_inst-1];
      n = pylab.fromfile(fileobj, 'i', 1)[0]
      taille = int(cplx*nbx);
      # we read values related to the required snapshot
      fileobj.seek(offset*prec_size, 1);
      data = pylab.fromfile(fileobj, prec, taille)
      if (cplx == 2):
        M = data[0:taille:2] + 1j*data[1:taille:2];
      else:
        M = data;
      
      fileobj.close()
      
      U = M;
      V = U;
      if (raff != 1):
        V = Refine1D(0.0, 1.0, U, raff)

      X = pylab.linspace(xmin, xmax, (nbx-1)*raff+1);
      Y = pylab.linspace(ymin, ymax, (nbx-1)*raff+1);
      XI = X; YI = Y;
      
    elif (type_grid == 2):
      # SismoPoint
      offset = 0;
      xmin_ = pylab.zeros(nb_grids);
      ymin_ = pylab.zeros(nb_grids);
      for i in range(nb_grids):
        xmin_[i] = pylab.fromfile(fileobj, prec, 1)[0]
        ymin_[i] = pylab.fromfile(fileobj, prec, 1)[0]
        if (i < num_inst-1):
          offset = offset + 1
        
      xmin = xmin_[num_inst-1];
      ymin = ymin_[num_inst-1];
      n = pylab.fromfile(fileobj, 'i', 1)[0]
      taille = cplx;
      # we read values related to the required snapshot
      fileobj.seek(offset*prec_size, 1);
      data = pylab.fromfile(fileobj, prec, taille)
      if (cplx == 2):
        M = data[0:taille:2] + 1j*data[1:taille:2];
      else:
        M = data;
      
      fileobj.close()
      
      XI = xmin
      YI = ymin
      V = M
      
    elif (type_grid == 3):
      
      # SismoCircle      
      offset = 0;
      xmin_ = pylab.zeros(nb_grids);
      ymin_ = pylab.zeros(nb_grids);
      radius_ = pylab.zeros(nb_grids);
      nbx_ = pylab.zeros(nb_grids).astype('int');
      for i in range(nb_grids):
        xmin_[i] = pylab.fromfile(fileobj, prec, 1)[0]
        ymin_[i] = pylab.fromfile(fileobj, prec, 1)[0]
        radius_[i] = pylab.fromfile(fileobj, prec, 1)[0]
        nbx_[i] = pylab.fromfile(fileobj, 'i', 1)[0]
        if (i < num_inst-1):
          offset = offset + nbx_[i];
        
      
      xmin = xmin_[num_inst-1];
      ymin = ymin_[num_inst-1];
      radius = radius_[num_inst-1];
      nbx = nbx_[num_inst-1];
      n = pylab.fromfile(fileobj, 'i', 1)[0]
      taille = int(cplx*nbx);
      # we read values related to the required snapshot
      fileobj.seek(offset*prec_size, 1);
      data = pylab.fromfile(fileobj, prec, taille)
      if (cplx == 2):
        M = data[0:taille:2] + 1j*data[1:taille:2];
      else:
        M = data;
      
      fileobj.close()
      
      teta = pylab.linspace(0.0, 2.0*pylab.pi, (nbx-1)*raff+1);
      X = radius*cos(teta) + xmin
      Y = radius*sin(teta) + ymin
      U = M
      V = U;
      if (raff != 1):
        V = Refine1D(0.0, 2.0*pylab.pi, U, raff)
      
      XI = X; YI = Y;
      
    elif (type_grid == 6):
      
      # SismoPointsFile
      offset = 0;
      nbx_ = pylab.zeros(nb_grids, dtype='int');
      for i in range(nb_grids):
        nbx_[i] = pylab.fromfile(fileobj, 'i', 1)[0]
        coord = pylab.fromfile(fileobj, prec, 2*int(nbx_[i]))
        if (i == num_inst-1):
          pts = coord
        
        if (i < num_inst-1):
          offset = offset + nbx_[i];          
        
      nbx = int(nbx_[num_inst-1]);
      n = pylab.fromfile(fileobj, 'i', 1)[0]
      taille = int(cplx*nbx);
      # we read values related to the required snapshot
      fileobj.seek(offset*prec_size, 1);
      data = pylab.fromfile(fileobj, prec, taille)
      if (cplx == 2):
        M = data[0:taille:2] + 1j*data[1:taille:2];
      else:
        M = data;
      
      fileobj.close()
      
      XI = pts[range(0,2*nbx,2)]
      YI = pts[range(1,2*nbx,2)]
      V = M
            
  else:
    if (type_grid == 0):
      # SismoGrille3D
      offset = 0;
      xmin_ = pylab.zeros(nb_grids); xmax_ = pylab.zeros(nb_grids)
      ymin_ = pylab.zeros(nb_grids); ymax_ = pylab.zeros(nb_grids)
      zmin_ = pylab.zeros(nb_grids); zmax_ = pylab.zeros(nb_grids)
      nbx_ = list(range(nb_grids)); nby_ = list(range(nb_grids)); nbz_ = list(range(nb_grids));
      for i in range(nb_grids):
        xmin_[i] = pylab.fromfile(fileobj, prec, 1)[0]
        xmax_[i] = pylab.fromfile(fileobj, prec, 1)[0]
        nbx_[i] = pylab.fromfile(fileobj, 'i', 1)[0]
        ymin_[i] = pylab.fromfile(fileobj, prec, 1)[0]
        ymax_[i] = pylab.fromfile(fileobj, prec, 1)[0]
        nby_[i] = pylab.fromfile(fileobj, 'i', 1)[0]
        zmin_[i] = pylab.fromfile(fileobj, prec, 1)[0]
        zmax_[i] = pylab.fromfile(fileobj, prec, 1)[0]
        nbz_[i] = pylab.fromfile(fileobj, 'i', 1)[0]
        if (i < num_inst-1):
          offset = offset + nbx_[i]*nby_[i]*nbz_[i];
      
      n = pylab.fromfile(fileobj, 'i', 1)[0]
      xmin = xmin_[num_inst-1]; xmax = xmax_[num_inst-1]; nbx = nbx_[num_inst-1];
      ymin = ymin_[num_inst-1]; ymax = ymax_[num_inst-1]; nby = nby_[num_inst-1];
      zmin = zmin_[num_inst-1]; zmax = zmax_[num_inst-1]; nbz = nbz_[num_inst-1];
      taille = cplx*(nbx*nby*nbz)
      # we read values related to the required snapshot
      fileobj.seek(offset*prec_size, 1);
      data = pylab.fromfile(fileobj, prec, taille)
      if (cplx == 2):
        M = data[0:taille:2] + 1j*data[1:taille:2];
      else:
        M = data;
      
      fileobj.close()
      
      X = pylab.linspace(xmin, xmax, nbx);
      Y = pylab.linspace(ymin, ymax, nby);
      Z = pylab.linspace(zmin, zmax, nbz);
      if (raff != 1):
        print("Refinement not implemented in 3-D")
      
      # here we obtain a 3-D array with inverted x and z
      V = pylab.reshape(M, [nbz, nby, nbx]);
      XI = X; YI = Y; ZI = Z;
      
    elif (type_grid == 1):
      # SismoGrille
      offset = 0;
      xmin_ = pylab.zeros(nb_grids); xmax_ = pylab.zeros(nb_grids);
      ymin_ = pylab.zeros(nb_grids); ymax_ = pylab.zeros(nb_grids);
      zmin_ = pylab.zeros(nb_grids); zmax_ = pylab.zeros(nb_grids);
      nbx_ = pylab.zeros(nb_grids, dtype='int');
      nby_ = pylab.zeros(nb_grids, dtype='int');
      nbz_ = pylab.zeros(nb_grids, dtype='int');
      coorx_ = pylab.zeros(nb_grids);
      coory_ = pylab.zeros(nb_grids);
      coorz_ = pylab.zeros(nb_grids);
      for i in range(nb_grids):
        xmin_[i] = pylab.fromfile(fileobj, prec, 1)[0]
        xmax_[i] = pylab.fromfile(fileobj, prec, 1)[0]
        nbx_[i] = pylab.fromfile(fileobj, 'i', 1)[0]
        ymin_[i] = pylab.fromfile(fileobj, prec, 1)[0]
        ymax_[i] = pylab.fromfile(fileobj, prec, 1)[0]
        nby_[i] = pylab.fromfile(fileobj, 'i', 1)[0]
        zmin_[i] = pylab.fromfile(fileobj, prec, 1)[0]
        zmax_[i] = pylab.fromfile(fileobj, prec, 1)[0]
        nbz_[i] = pylab.fromfile(fileobj, 'i', 1)[0]
        coorx_[i] = pylab.fromfile(fileobj, prec, 1)[0]
        coory_[i] = pylab.fromfile(fileobj, prec, 1)[0]
        coorz_[i] = pylab.fromfile(fileobj, prec, 1)[0]
        if (i < num_inst-1):
          offset = offset + nbx_[i]*nby_[i] + nbx_[i]*nbz_[i] + nby_[i]*nbz_[i];
      
      xmin = xmin_[num_inst-1]; xmax = xmax_[num_inst-1];
      ymin = ymin_[num_inst-1]; ymax = ymax_[num_inst-1];
      zmin = zmin_[num_inst-1]; zmax = zmax_[num_inst-1];
      coor = pylab.zeros(3); coor[0] = coorx_[num_inst-1];
      coor[1] = coory_[num_inst-1]; coor[2] = coorz_[num_inst-1];
      nbx = nbx_[num_inst-1]; nby = nby_[num_inst-1]; nbz = nbz_[num_inst-1];
      n = pylab.fromfile(fileobj, 'i', 1)[0]
      taille = int(cplx*(nbx*nby + nbx*nbz + nby*nbz));
      # we read values related to the required snapshot
      fileobj.seek(offset*prec_size, 1);
      data = pylab.fromfile(fileobj, prec, taille)
      if (cplx == 2):
        M = data[0:taille:2] + 1j*data[1:taille:2];
      else:
        M = data;
      
      fileobj.close()

      V1 = pylab.reshape(M[0:nby*nbz], (nby, nbz));
      V2 = pylab.reshape(M[nby*nbz:nbx*nbz+nby*nbz], (nbx, nbz));
      V3 = pylab.reshape(M[nbz*(nbx+nby):nbz*(nbx+nby)+nbx*nby], (nbx, nby));
      if (raff != 1):
        V1 = Refine2D(ymin, ymax, zmin, zmax, V1, raff)
        V2 = Refine2D(xmin, xmax, zmin, zmax, V2, raff)
        V3 = Refine2D(xmin, xmax, ymin, ymax, V3, raff)

      X = pylab.linspace(xmin, xmax, (nbx-1)*raff+1);
      Y = pylab.linspace(ymin, ymax, (nby-1)*raff+1);
      Z = pylab.linspace(zmin, zmax, (nbz-1)*raff+1);
      XI = X
      YI = Y
      ZI = Z
      V = [pylab.transpose(V1),pylab.transpose(V2),pylab.transpose(V3)];
    elif (type_grid == 2):
      # SismoPlane
      offset = 0;
      xmin_ = pylab.zeros(nb_grids); xmax_ = pylab.zeros(nb_grids);
      ymin_ = pylab.zeros(nb_grids); ymax_ = pylab.zeros(nb_grids);
      zmin_ = pylab.zeros(nb_grids); zmax_ = pylab.zeros(nb_grids);
      centerx_ = pylab.zeros(nb_grids); centery_ = pylab.zeros(nb_grids);
      centerz_ = pylab.zeros(nb_grids);
      nbx_ = pylab.zeros(nb_grids, dtype='int'); nby_ = pylab.zeros(nb_grids, dtype='int');
      for i in range(nb_grids):
        xmin_[i] = pylab.fromfile(fileobj, prec, 1)[0]
        xmax_[i] = pylab.fromfile(fileobj, prec, 1)[0]
        ymin_[i] = pylab.fromfile(fileobj, prec, 1)[0]
        ymax_[i] = pylab.fromfile(fileobj, prec, 1)[0]
        zmin_[i] = pylab.fromfile(fileobj, prec, 1)[0]
        zmax_[i] = pylab.fromfile(fileobj, prec, 1)[0]
        centerx_[i] = pylab.fromfile(fileobj, prec, 1)[0]
        centery_[i] = pylab.fromfile(fileobj, prec, 1)[0]
        centerz_[i] = pylab.fromfile(fileobj, prec, 1)[0]
        nbx_[i] = pylab.fromfile(fileobj, 'i', 1)[0]
        nby_[i] = pylab.fromfile(fileobj, 'i', 1)[0]
        if (i < num_inst-1):
          offset = offset + nbx_[i]*nby_[i];
      
      origin = pylab.zeros(3); extA = pylab.zeros(3); extB = pylab.zeros(3)
      origin[0] = xmin_[num_inst-1]; extA[0] = xmax_[num_inst-1];
      origin[1] = ymin_[num_inst-1]; extA[1] = ymax_[num_inst-1];
      origin[2] = zmin_[num_inst-1]; extA[2] = zmax_[num_inst-1];
      extB[0] = centerx_[num_inst-1];
      extB[1] = centery_[num_inst-1]; extB[2] = centerz_[num_inst-1];
      nbx = nbx_[num_inst-1]; nby = nby_[num_inst-1];
      n = pylab.fromfile(fileobj, 'i', 1)[0]
      taille = cplx*nbx*nby;
      # we read values related to the required snapshot
      fileobj.seek(offset*prec_size, 1);
      data = pylab.fromfile(fileobj, prec, taille)
      if (cplx == 2):
        M = data[0:taille:2] + 1j*data[1:taille:2];
      else:
        M = data;
      
      fileobj.close()

      # we form the transpose of the original field (for display purposes)
      Tx = pylab.linspace(0, 1, (nbx-1)*raff+1);
      Ty = pylab.linspace(0, 1, (nby-1)*raff+1);
      [Tx, Ty] = pylab.meshgrid(Tx, Ty);
      U = pylab.reshape(M, (nby, nbx));
      if (raff != 1):
        U = Refine2D(0.0, 1.0, 0.0, 1.0, U, raff)
      
      V = U
      
      XI = (extA[0]-origin[0])*Tx + (extB[0]-origin[0])*Ty + origin[0];
      YI = (extA[1]-origin[1])*Tx + (extB[1]-origin[1])*Ty + origin[1];
      ZI = (extA[2]-origin[2])*Tx + (extB[2]-origin[2])*Ty + origin[2];
      
    elif (type_grid == 3):
      # SismoLine
      offset = 0;
      xmin_ = pylab.zeros(nb_grids); xmax_ = pylab.zeros(nb_grids);
      ymin_ = pylab.zeros(nb_grids); ymax_ = pylab.zeros(nb_grids);
      zmin_ = pylab.zeros(nb_grids); zmax_ = pylab.zeros(nb_grids);
      nbx_ = list(range(nb_grids))
      for i in range(nb_grids):
        xmin_[i] = pylab.fromfile(fileobj, prec, 1)[0]
        xmax_[i] = pylab.fromfile(fileobj, prec, 1)[0]
        ymin_[i] = pylab.fromfile(fileobj, prec, 1)[0]
        ymax_[i] = pylab.fromfile(fileobj, prec, 1)[0]
        zmin_[i] = pylab.fromfile(fileobj, prec, 1)[0]
        zmax_[i] = pylab.fromfile(fileobj, prec, 1)[0]
        nbx_[i] = pylab.fromfile(fileobj, 'i', 1)[0]
        if (i < num_inst-1):
          offset = offset + nbx_[i];
      
      xmin = xmin_[num_inst-1]; xmax = xmax_[num_inst-1];
      ymin = ymin_[num_inst-1]; ymax = ymax_[num_inst-1];
      zmin = zmin_[num_inst-1]; zmax = zmax_[num_inst-1];
      nbx = nbx_[num_inst-1];
      n = pylab.fromfile(fileobj, 'i', 1)[0]
      taille = cplx*nbx;
      # we read values related to the required snapshot
      fileobj.seek(offset*prec_size, 1);
      data = pylab.fromfile(fileobj, prec, taille)
      if (cplx == 2):
        M = data[0:taille:2] + 1j*data[1:taille:2];
      else:
        M = data;
      
      fileobj.close()
      
      T = pylab.linspace(0, 1, (nbx-1)*raff + 1); 
      V = M
      if (raff != 1):
        V = Refine1D(0.0, 1.0, M, raff)

      XI = pylab.linspace(xmin, xmax, nbx);
      YI = pylab.linspace(ymin, ymax, nbx);
      ZI = pylab.linspace(zmin, zmax, nbx);
      
    elif (type_grid == 4):
      # SismoPoint
      offset = 0;
      xmin_ = pylab.zeros(nb_grids);
      ymin_ = pylab.zeros(nb_grids);
      zmin_ = pylab.zeros(nb_grids);
      for i in range(nb_grids):
        xmin_[i] = pylab.fromfile(fileobj, prec, 1)[0]
        ymin_[i] = pylab.fromfile(fileobj, prec, 1)[0]
        zmin_[i] = pylab.fromfile(fileobj, prec, 1)[0]
        if (i < num_inst-1):
          offset = offset + 1;
      
      xmin = xmin_[num_inst-1];
      ymin = ymin_[num_inst-1];
      zmin = zmin_[num_inst-1];
      n = pylab.fromfile(fileobj, 'i', 1)[0]
      taille = cplx;
      # we read values related to the required snapshot
      fileobj.seek(offset*prec_size, 1);
      data = pylab.fromfile(fileobj, prec, taille)
      if (cplx == 2):
        M = data[0:taille:2] + 1j*data[1:taille:2];
      else:
        M = data;
      
      fileobj.close()
      
      XI = xmin;
      YI = ymin;
      ZI = zmin;
      V = M;
      
    elif (type_grid == 5):
      # SismoPointsFile
      offset = 0;
      nbx_ = list(range(nb_grids))
      for i in range(nb_grids):
        nbx_[i] = pylab.fromfile(fileobj, 'i', 1)[0]
        coord = pylab.fromfile(fileobj, prec, 3*nbx_[i])
        if (i == num_inst-1):
          pts = coord
        if (i < num_inst-1):
          offset = offset + nbx_[i];

      nbx = int(nbx_[num_inst-1]);
      n = pylab.fromfile(fileobj, 'i', 1)[0]
      taille = cplx*nbx;
      # we read values related to the required snapshot
      fileobj.seek(offset*prec_size, 1);
      data = pylab.fromfile(fileobj, prec, taille)
      if (cplx == 2):
        M = data[0:taille:2] + 1j*data[1:taille:2];
      else:
        M = data;
      
        fileobj.close()
      
      XI = pts[range(0,3*nbx,3)]
      YI = pts[range(1,3*nbx,3)]
      ZI = pts[range(2,3*nbx,3)]
      V = M
    elif (type_grid == 6):
      
      # SismoCircle
      offset = 0;
      centerx_ = pylab.zeros(nb_grids)
      centery_ = pylab.zeros(nb_grids)
      centerz_ = pylab.zeros(nb_grids)
      normalex_ = pylab.zeros(nb_grids)
      normaley_ = pylab.zeros(nb_grids)
      normalez_ = pylab.zeros(nb_grids)
      radiusA_ = pylab.zeros(nb_grids)
      radiusB_ = pylab.zeros(nb_grids)
      nbx_ = list(range(nb_grids))
      for i in range(nb_grids):
        centerx_[i] = pylab.fromfile(fileobj, prec, 1)[0]
        centery_[i] = pylab.fromfile(fileobj, prec, 1)[0]
        centerz_[i] = pylab.fromfile(fileobj, prec, 1)[0]
        normalex_[i] = pylab.fromfile(fileobj, prec, 1)[0]
        normaley_[i] = pylab.fromfile(fileobj, prec, 1)[0]
        normalez_[i] = pylab.fromfile(fileobj, prec, 1)[0]
        radiusA_[i] = pylab.fromfile(fileobj, prec, 1)[0]
        radiusB_[i] = pylab.fromfile(fileobj, prec, 1)[0]
        nbx_[i] = pylab.fromfile(fileobj, 'i', 1)[0]
        if (i < num_inst-1):
          offset = offset + nbx_[i];

      centerx = centerx_[num_inst-1]; radiusA = radiusA_[num_inst-1];
      centery = centery_[num_inst-1]; radiusB = radiusB_[num_inst-1];
      centerz = centerz_[num_inst-1];
      normale = pylab.array([normalex_[num_inst-1], normaley_[num_inst-1], normalez_[num_inst-1]])
      nbx = nbx_[num_inst-1];
      n = pylab.fromfile(fileobj, 'i', 1)[0]
      taille = cplx*nbx;
      # we read values related to the required snapshot
      fileobj.seek(offset*prec_size, 1);
      data = pylab.fromfile(fileobj, prec, taille)
      if (cplx == 2):
        M = data[0:taille:2] + 1j*data[1:taille:2];
      else:
        M = data;
      
      fileobj.close()

      # finding the two vectors of plane, knowing the normale
      perm = pylab.array([0, 1, 2])
      if (abs(normale[0]) > max(abs(normale[1]),abs(normale[2]))):
        perm[2] = 0
        if (abs(normale[1]) > abs(normale[2])):
          perm[1] = 1
          perm[0] = 2
        else:
          perm[1] = 2
          perm[0] = 1
      else:
        if (abs(normale[1]) > abs(normale[2])):
          perm[2] = 1
          if (abs(normale[2]) > abs(normale[0])):
            perm[1] = 2
            perm[0] = 0
          else:
            perm[1] = 0
            perm[0] = 2
        else:
          perm[2] = 2
          if (abs(normale[1]) > abs(normale[0])):
            perm[1] = 1
            perm[0] = 0
          else:
            perm[1] = 0
            perm[0] = 1
      
      vec_u = pylab.zeros(3); vec_u[perm[2]] = -normale[perm[1]];
      vec_u[perm[1]] = normale[perm[2]]
      vec_v = pylab.cross(normale, vec_u)
      
      V = M
      if (raff != 1):
        V = Refine1D(0.0, 1.0, M, raff)

      teta = pylab.linspace(0, 2*pi, (nbx-1)*raff+2); teta = teta[0:(nbx-1)*raff+1];
      XI = radiusA*vec_u[0]*cos(teta) + radiusB*vec_v[0]*sin(teta) + centerx;
      YI = radiusA*vec_u[1]*cos(teta) + radiusB*vec_v[1]*sin(teta) + centery;
      ZI = radiusA*vec_u[2]*cos(teta) + radiusB*vec_v[2]*sin(teta) + centerz;
      
        
  return XI,YI,ZI,coor,V

def loadEH(nom_fichier, num_inst=1):
  """ Reading of a "Montjoie" output file (MATLAB output) 
  usage : [X, Y, Z, coor, Ex, Ey, Ez, Hx, Hy, Hz] = loadEH(nom_fichier)
          [X, Y, Z, coor, Ex, Ey, Ez, Hx, Hy, Hz] = loadEH(nom_fichier, num_inst) """
  [X, Y, Z, coor, Ex] = loadND(nom_fichier + "_U0.dat")
  [X, Y, Z, coor, Ey] = loadND(nom_fichier + "_U1.dat")
  [X, Y, Z, coor, Ez] = loadND(nom_fichier + "_U2.dat")
  [X, Y, Z, coor, Hx] = loadND(nom_fichier + "_dU0.dat")
  [X, Y, Z, coor, Hy] = loadND(nom_fichier + "_dU1.dat")
  [X, Y, Z, coor, Hz] = loadND(nom_fichier + "_dU2.dat")
  return X, Y, Z, coor, Ex, Ey, Ez, Hx, Hy, Hz

def loadE(nom_fichier, num_inst=1):
  """ Reading of a "Montjoie" output file (MATLAB output) 
  usage : [X, Y, Z, coor, Ex, Ey, Ez, Hx, Hy, Hz] = loadEH(nom_fichier)
          [X, Y, Z, coor, Ex, Ey, Ez, Hx, Hy, Hz] = loadEH(nom_fichier, num_inst) """
  [X, Y, Z, coor, Ex] = loadND(nom_fichier + "_U0.dat")
  [X, Y, Z, coor, Ey] = loadND(nom_fichier + "_U1.dat")
  [X, Y, Z, coor, Ez] = loadND(nom_fichier + "_U2.dat")
  return X, Y, Z, coor, Ex, Ey, Ez

def film2D(root, ext, n1, n2, cmin, cmax, num_inst = 1, raff=1, dt=0, display_nan=True):
  """ Displays 2-D solutions, which are stored in files 
     of the type toto0000Ext.dat, toto0001Ext.dat, toto0002Ext.dat 

  Usage : film2D('toto', 'Ext', 0, 500, -0.01, 0.01) 
          film2D(base, ext, n1, n2, cmin, cmax) """
  
  h = pylab.figure()
  pylab.ion()
  pylab.draw()
  nom = root + EntierToString(n1) + ext
  [X, Y, Z, coor, V] = loadND(nom, num_inst, raff)
  if (display_nan):
    num, = pylab.nonzero(pylab.ravel(abs(V)==0))
    n = V.shape[1]
    V[num//n, num%n] = pylab.nan

  xmin = min(X); xmax = max(X); ymin = min(Y); ymax = max(Y)
  hh = pylab.imshow(V, vmin=cmin, vmax=cmax, interpolation='bilinear', origin='lower', extent=(xmin, xmax, ymin, ymax))
  if (dt != 0):
    os.system('sleep ' +str(dt))
  for i in range(n1+1, n2):
    nom = root + EntierToString(i) + ext
    [X, Y, Z, coor, V] = loadND(nom, num_inst, raff)
    if (display_nan):
      num, = pylab.nonzero(pylab.ravel(abs(V)==0))
      n = V.shape[1]
      V[num//n, num%n] = pylab.nan
    
    hh.set_data(V)
    pylab.draw()
    if (dt != 0):
      os.system('sleep ' +str(dt))

def loadAll(racine, ext, debut, fin, type_entier = 0, taille = 0):
  """ Loads a list of files stored in loadND format (1-D solutions only)
  usage : t, V = loadAll(racine, extension, debut, fin)
          t, V = loadAll(racine, extension, debut, fin, type_entier)
   
   racine : root of the output files
   ext : extension of the output files
   debut, fin : files from racine debut ext.dat until racine fin ext.dat are loaded
   type_entier : if equal to 0, the numbers are written in 4 characters
   returns t (X parameter returned by loadND, assuming that all the files give the same X)
   returns V (V parameter returned by loadND, stored as a matrix)
   
   for example if racine = 'test', ext = '.dat', debut =0, and fin = 10, the files
   test0000.dat, test0001.dat, ..., test0009.dat will be loaded (if type_entier = 0) """

  if (type_entier == 0):
    [X, Y, Z, coor, V] = loadND(racine+EntierToString(fin-1)+ext)
  else:
    [X, Y, Z, coor, V] = loadND(racine+str(fin-1)+ext)
    
  t = X
  if (taille == 0):
    taille = V.size
  else:
    dt = t[1]-t[0]
    t = pylab.linspace(t[0], (taille-1)*dt + t[0], taille)
  
  E = pylab.zeros([taille, fin-debut]) + 1j*pylab.zeros([taille, fin-debut])    
  for i in range(debut, fin):
    if (type_entier == 0):
      [X, Y, Z, coor, V] = loadND(racine+EntierToString(i)+ext)
    else:
      [X, Y, Z, coor, V] = loadND(racine+str(i)+ext)
      
    test = True
    if (V.size == E.shape[0]):
      n0 = 0
      n1 = V.size
    elif (V.size < E.shape[0]):
      n0 = (E.shape[0]-V.size)/2
      n1 = n0 + V.size
    else:
      print("Error : the data is too large ")
      print("Tailles", V.size, E.shape[0])
      test = False
    
    if (test):
      E[n0:n1, i-debut] = V

  return t, E


def GetLocalMaxima(V, threshold = 0):
  """ returns the position of local maxima in vector V
  usage : num = GetLocalMaxima(V) """
  N = V.size
  diff = abs(V[0:N-1]) - abs(V[1:N])
  prod = diff[0:N-2]*diff[1:N-1]
  num = pylab.nonzero(pylab.ravel(prod < 0))[0] + 1
  n = []
  for p in num:
    if abs(V[p]) > threshold:
      n.append(p)
      
  return n

def GetEnvelope(V, N0 = 0, N1 = 0):
  """ returns the envelope of a signal V
    usage : Venv = GetEnvelope(V)
            Venv = GetEnvelope(V, Ntronc)
    V : data (it can be a vector or a matrix)
    N0, N1 : only fft between N0 and N1 is kept
    If V is a matrix, the function is applied to each column of V """
  
  threshold = 1e-6
  if (len(V.shape) == 1):
    N = V.shape[0]
    # Fourier transform is computed
    Vchap = pylab.fft(V)

    # we keep only the signal that will contain the envelope
    if ((N0 == 0) and (N1 == 0)):
      num = GetLocalMaxima(Vchap, threshold*abs(Vchap).max())
      deb = (num[0] + num[1])//2
      Vchap[deb:] = 0
    else:
      Vchap[0:N0] = 0
      Vchap[N1:] = 0
    
    # returning back to real space
    Venv = pylab.ifft(Vchap)
    return Venv
  elif (len(V.shape) == 2):
    N = V.shape[0]
    Venv = pylab.zeros(V.shape) + 1j*pylab.zeros(V.shape)
    for i in range(V.shape[1]):
      # Fourier transform is computed
      Vchap = pylab.fft(V[:,i])

      # we keep only the signal that will contain the envelope
      if ((N0 == 0) and (N1 == 0)):
        num = GetLocalMaxima(Vchap, threshold*abs(Vchap).max())
        deb = (num[0] + num[1])//2
        Vchap[deb:] = 0
      else:
        Vchap[0:N0] = 0
        Vchap[N1:] = 0

      # returning back to real space
      Venv[:,i] = pylab.ifft(Vchap)
    
    return Venv
  else:
    print("case not handled")

    
def plot2dinst(X, Y, V, cmin='auto', cmax='auto', display_nan = True, log_scale = False,
               display_colorbar = True, aspect_fig=None, cmap = None, hide_axis=False):
  """ Displays a 2-D solution by using imshow
  usage : plot2dinst(X, Y, V)
          plot2dinst(X, Y, real(V))
          plot2dinst(X, Y, V, cmin, cmax)
          If you have a matrix, you can generate X, Y with the following command
          X = pylab.linspace(xmin, xmax, nbx)
          Y = pylab.linspace(ymin, ymax, nby)
          
  X : x-coordinates as given by function loadND
  Y : y-coordinates as given by function loadND
  V : values of the 2-D solution
  cmin , cmax : the colorbar is set for values between cmin and cmax """
  
  xmin = min(pylab.reshape(X,(pylab.size(X))))
  ymin = min(pylab.reshape(Y,(pylab.size(Y))))
  xmax = max(pylab.reshape(X,(pylab.size(X))))
  ymax = max(pylab.reshape(Y,(pylab.size(Y))))
  if (cmin == 'auto'):
    y = pylab.reshape(V, pylab.size(V))
    cmin = min(y)
    cmax = max(y)
  
  V2 = V.copy()
  if (display_nan):
    num, = pylab.nonzero(pylab.ravel(abs(V)==0))
    if (len(num) > 0):
      n = V.shape[1]
      V2[num//n, num%n] = pylab.nan

  color_map = cmap
  if (cmap == None):
    color_map = pylab.cm.jet
  
  h = pylab.figure();
  if (aspect_fig == None):
    if (log_scale):
      hh = pylab.imshow(V2, interpolation='bilinear', origin='lower', extent=(xmin,xmax,ymin,ymax), norm=matplotlib.colors.LogNorm(vmin=cmin, vmax=cmax), cmap=color_map)
    else:
      hh = pylab.imshow(V2, vmin=cmin, vmax=cmax, interpolation='bilinear', origin='lower', extent=(xmin,xmax,ymin,ymax), cmap=color_map)
    pylab.axis([xmin, xmax, ymin, ymax])    
    pylab.axis('image')
    if (hide_axis):
      pylab.axis('off');
      h.axes[0].get_xaxis().set_visible(False)
      h.axes[0].get_yaxis().set_visible(False)
  else:
    hh = pylab.imshow(V2, vmin=cmin, vmax=cmax, interpolation='bilinear', origin='lower', extent=(xmin,xmax,ymin,ymax), aspect='auto', cmap=color_map)

  if (display_colorbar):
    if (log_scale):
      #tab = 10.0**pylab.arange(pylab.ceil(pylab.log10(cmin)), pylab.floor(pylab.log10(cmax))+1.0, 0.5);
      #pylab.colorbar(format="%.3g", ticks = tab)
      pylab.colorbar(format="%.3g")
    else:
      pylab.colorbar(format="%.3g")
  
  return hh


def plot3d_plane(X, Y, Z, V, coor, cmin = 'auto', cmax = 'auto', component = 0, display_nan = True, cmap = None):
  """ Displays a 3-D solution by using imshow for the three planes
  usage : plot3d_plane(X, Y, Z, V, coor)
          plot3d_plane(X, Y, Z, real(V))
          plot3d_plane(X, Y, Z, V, coor, cmin, cmax)
          
  X : x-coordinates as given by function loadND
  Y : y-coordinates as given by function loadND
  Z : z-coordinates as given by function loadND
  coor : intersection of three planes as given by function loadND
  cmin , cmax : the colorbar is set for values between cmin and cmax """
  
  xmin = min(pylab.reshape(X,(pylab.size(X))))
  ymin = min(pylab.reshape(Y,(pylab.size(Y))))
  zmin = min(pylab.reshape(Z,(pylab.size(Z))))
  xmax = max(pylab.reshape(X,(pylab.size(X))))
  ymax = max(pylab.reshape(Y,(pylab.size(Y))))
  zmax = max(pylab.reshape(Z,(pylab.size(Z))))
  h = pylab.figure()  
  
  pylab.subplot(2,2,2)
  V0 = V[0].copy()
  V1 = V[1].copy()
  V2 = V[2].copy()
  if (cmin=='auto'):
    cmin = min(pylab.reshape(V0,(pylab.size(V0),1)))[0]
    cmax = max(pylab.reshape(V0,(pylab.size(V0),1)))[0]
    cmin = min(cmin, min(pylab.reshape(V1,(pylab.size(V1),1)))[0])
    cmax = max(cmax, max(pylab.reshape(V1,(pylab.size(V1),1)))[0])
    cmin = min(cmin, min(pylab.reshape(V2,(pylab.size(V2),1)))[0])
    cmax = max(cmax, max(pylab.reshape(V2,(pylab.size(V2),1)))[0])

  color_map = cmap
  if (cmap == None):
    color_map = pylab.cm.jet
  
  if (display_nan):
    num, = pylab.nonzero(pylab.ravel(abs(V0)==0))
    n = V0.shape[1]
    V0[num//n, num%n] = pylab.nan
    num, = pylab.nonzero(pylab.ravel(abs(V1)==0))
    n = V1.shape[1]
    V1[num//n, num%n] = pylab.nan
    num, = pylab.nonzero(pylab.ravel(abs(V2)==0))
    n = V2.shape[1]
    V2[num//n, num%n] = pylab.nan

  if (component == 0):
    pylab.imshow(pylab.real(V0), vmin=cmin, vmax=cmax, interpolation='bilinear', origin='lower', extent=(ymin,ymax,zmin,zmax), cmap=color_map)
    pylab.axis('equal')
    
    pylab.subplot(2,2,3)
    pylab.imshow(pylab.real(V1), vmin=cmin, vmax=cmax, interpolation='bilinear', origin='lower', extent=(xmin,xmax,zmin,zmax), cmap=color_map)
    pylab.axis('equal')
    
    pylab.subplot(2,2,4)
    pylab.imshow(pylab.real(V2), vmin=cmin, vmax=cmax, interpolation='bilinear', origin='lower', extent=(xmin,xmax,ymin,ymax), cmap=color_map)
    pylab.axis('equal')
  elif (component == 1):
    pylab.imshow(pylab.imag(V0), vmin=cmin, vmax=cmax, interpolation='bilinear', origin='lower', extent=(ymin,ymax,zmin,zmax), cmap=color_map)
    pylab.axis('equal')
    
    pylab.subplot(2,2,3)
    pylab.imshow(pylab.imag(V1), vmin=cmin, vmax=cmax, interpolation='bilinear', origin='lower', extent=(xmin,xmax,zmin,zmax), cmap=color_map)
    pylab.axis('equal')
    
    pylab.subplot(2,2,4)
    pylab.imshow(pylab.imag(V2), vmin=cmin, vmax=cmax, interpolation='bilinear', origin='lower', extent=(xmin,xmax,ymin,ymax), cmap=color_map)
    pylab.axis('equal')
  else:
    pylab.imshow(abs(V0), vmin=cmin, vmax=cmax, interpolation='bilinear', origin='lower', extent=(ymin,ymax,zmin,zmax), cmap=color_map)
    pylab.axis('equal')
    
    pylab.subplot(2,2,3)
    pylab.imshow(abs(V1), vmin=cmin, vmax=cmax, interpolation='bilinear', origin='lower', extent=(xmin,xmax,zmin,zmax), cmap=color_map)
    pylab.axis('equal')
    
    pylab.subplot(2,2,4)
    pylab.imshow(abs(V2), vmin=cmin, vmax=cmax, interpolation='bilinear', origin='lower', extent=(xmin,xmax,ymin,ymax), cmap=color_map)
    pylab.axis('equal')

  pylab.subplot(2,2,1)
  pylab.axis('off')
  pylab.colorbar(format="%.3g")

def erreurL2(*args, **kwargs):
  """ Returns L2 error between two numeric arrays V and W
  Usage : err = erreurL2(V, W)
          err = erreurL2(V, W, threshold=epsilon)
          err = erreurL2(Ux, UxRef, Uy, UyRef, ...)
  
  V : first numeric array (it can be either a vector, matrix or 3-D array)
  W : second numeric array (it can be either a vector, matrix or 3-D array)
  err : L2 error between the two vectors
  threshold : numerical values below threshold are not compared
  only simultaneous non-null values of V and W are considered """
  
  threshold = 1e-30
  for k, v in kwargs.items():
    if (k == 'threshold'):
      threshold = v
      
  err = 0.0
  norm_sol = 0.0
  for i in range(len(args)//2):
    x = pylab.reshape(args[2*i], pylab.size(args[2*i]))
    y = pylab.reshape(args[2*i+1], pylab.size(args[2*i+1]))
    droptol = threshold*max(abs(x))
    x = x*(abs(x) > droptol)*(abs(y) > droptol);
    y = y*(abs(x) > droptol)*(abs(y) > droptol);
    err += pylab.norm(x-y)**2
    norm_sol += pylab.norm(y)**2
    
  return pylab.sqrt(err/norm_sol)

def erreurMediane(V, W, threshold = 1e-30):
  """ Returns median error between two numeric arrays V and W
  Usage : err = erreurMediane(V, W)
          err = erreurMediane(V, W, threshold)
  
  V : first numeric array (it can be either a vector, matrix or 3-D array)
  W : second numeric array (it can be either a vector, matrix or 3-D array)
  err : median error between the two vectors
  threshold : numerical values below threshold are not compared
  only simultaneous non-null values of V and W are considered """
  
  ratio = 0.9
  x = pylab.reshape(V,(pylab.size(V)))
  y = pylab.reshape(W,(pylab.size(V)))
  droptol = threshold*max(abs(x))
  x = x*(abs(x) > droptol)*(abs(y) > droptol);
  y = y*(abs(x) > droptol)*(abs(y) > droptol);
  diff = pylab.sort(abs(x-y))
  val = pylab.sort(abs(y))
  offset = min(pylab.nonzero(pylab.ravel(val > 0))[0])
  N = pylab.size(x) - offset
  num = int(offset + ratio*N)
  err = diff[num]/val[num]
  return err

def erreurMode(V, mode):
  n=pylab.nonzero(pylab.ravel(abs(mode)== abs(mode).max()))[0][0]
  n1=n//V.shape[0]
  n2=n%V.shape[0]
  error = erreurL2(mode/mode[n1, n2], V/V[n1, n2])
  return error
        
def saveData(nom_fichier, V):
  """ Save an object into a file 
  Usage : saveData(nom_fichier, V) """
  
  f = open(nom_fichier,'w');
  pickle.dump(V, f);
  f.close();

def loadData(nom_fichier):
  """ Read an object from a file 
  Usage : V = loadData(nom_fichier) """  
  f = open(nom_fichier,'r');
  V = pickle.load(f);
  f.close();
  return V

def GetUniqueIndex(sig):
  """ Returns the indices index such that sig[index] is strictly increasing
  Usage : index = GetUniqueIndex(sig)
  sig : vector containing an almost increasing sequence (usually time ticks) """
  
  n = len(sig)
  Delta = sig[1:n] - sig[0:n-1]
  index, = pylab.nonzero(pylab.ravel(Delta<=0))
  if(len(index) == 0):
    b = range(n)
    return list(b)

  b = range(index[0]+1)
  for i in range(len(index)):
    iend, = pylab.nonzero(pylab.ravel(sig[index[i]+1:] > sig[index[i]-1]))
    fin = n-1
    if (i < len(index)-1):
      fin = index[i+1]
      
    if (len(iend) > 0):
      b2 = range(len(b) + fin - index[i]-iend[0]-1)
      b2[0:len(b)] = b
      b2[len(b):] = range(index[i]+iend[0]+2,fin+1)
      b = b2
    
  return b

def loadSismo(nom_fichier):
  """ loads a seismogramm and avoids reprises
  Usage : V = loadSismo(nom_fichier)
  The file is assumed to contain in the first column the times
  and values in the other columns. If the time column is not
  strictly increasing, increasing times are extracted with the help of GetUniqueIndex """
  V = pylab.loadtxt(nom_fichier)
  b = GetUniqueIndex(V[:,0])
  U = V[b,:]
  return U

def loadReverse(nom_fichier):
  V = pylab.loadtxt(nom_fichier)
  i, = pylab.nonzero(pylab.ravel(V[:,0] == V[:,0].min()))
  V[0:i+1, :] = pylab.flipud(V[0:i+1, :])
  return V
  
def load1D(nom_fichier):
  """ loads a complex vector (as produced by using Write function of Seldon vectors)
  Usage : V = load1D(nom_fichier)"""
  fileobj = open(nom_fichier, mode='rb')
  taille = np.fromfile(fileobj, 'i', 1)[0]
  Vz = np.fromfile(fileobj, 'd', 2*taille)
  V = Vz[0::2] + Vz[1::2]*1j;
  fileobj.close()
  return V

def load1D_real(nom_fichier):
  """ loads a real vector (as produced by using Write function of Seldon vectors)
  Usage : V = load1D_real(nom_fichier)"""
  fileobj = open(nom_fichier, mode='rb')
  taille = pylab.fromfile(fileobj, 'i', 1)[0]
  V = pylab.fromfile(fileobj, 'd', taille)
  fileobj.close()
  return V

def load1D_int(nom_fichier):
  """ loads a vector of integers (as produced by using Write function of Seldon vectors)
  Usage : V = load1D_int(nom_fichier)"""
  fileobj = open(nom_fichier, mode='rb')
  taille = pylab.fromfile(fileobj, 'i', 1)[0]
  V = pylab.fromfile(fileobj, 'i', taille)
  fileobj.close()
  return V

def load_full(nom_fichier, double_prec = True):
  """ loads a complex matrix (as produced by using Write function of Seldon dense matrices)
  Usage : A = load_full(nom_fichier)"""
  fileobj = open(nom_fichier, mode='rb')
  m = pylab.fromfile(fileobj, 'i', 1)[0]
  n = pylab.fromfile(fileobj, 'i', 1)[0]
  taille = m*n
  if (double_prec):
    Vz = pylab.fromfile(fileobj, 'd', 2*taille)
  else:
    Vz = pylab.fromfile(fileobj, 'f', 2*taille).astype('float64')
  
  V = Vz[0::2] + Vz[1::2]*1j;
  fileobj.close()
  A = pylab.reshape(V,(m, n))
  return A

def load_fullReal(nom_fichier, double_prec = True):
  """ loads a real matrix (as produced by using Write function of Seldon dense matrices)
  Usage : A = load_full(nom_fichier)"""
  fileobj = open(nom_fichier, mode='rb')
  m = pylab.fromfile(fileobj, 'i', 1)[0]
  n = pylab.fromfile(fileobj, 'i', 1)[0]
  taille = m*n
  if (double_prec):
    V = pylab.fromfile(fileobj, 'd', taille)
  else:
    V = pylab.fromfile(fileobj, 'f', taille).astype('float64')
  fileobj.close()
  A = pylab.reshape(V,(m, n))
  return A

def write_full(nom_fichier, A, double_prec = True):
  """ writes a real or complex matrix (so that it is readable by Seldon) 
  Usage : write_full(nom_fichier, A)
          write_full(nom_fichier, A, False) """
  fileobj = open(nom_fichier, mode='wb')
  taille = np.array(A.shape).astype('int32')
  taille.tofile(fileobj)
  if (double_prec):
    if ((A.dtype == 'float64') or (A.dtype == 'float32')):
      A.astype('float64').tofile(fileobj)
    elif ((A.dtype == 'int32') or (A.dtype == 'int64')):
      A.astype('int32').tofile(fileobj)
    else:
      A.astype('complex128').tofile(fileobj)
  else:
    if ((A.dtype == 'float64') or (A.dtype == 'float32')):
      A.astype('float32').tofile(fileobj)
    elif ((A.dtype == 'int32') or (A.dtype == 'int64')):
      A.astype('int32').tofile(fileobj)
    else:
      A.astype('complex64').tofile(fileobj)  

def load_fullSym(nom_fichier):
  """ loads a complex matrix (as produced by using Write function of Seldon dense symmetric matrices)
  Usage : A = load_full(nom_fichier)"""
  fileobj = open(nom_fichier, mode='rb')
  n = pylab.fromfile(fileobj, 'i', 1)[0]
  n = pylab.fromfile(fileobj, 'i', 1)[0]
  taille = (n+1)*n
  Vz = pylab.fromfile(fileobj, 'd', taille)
  V = Vz[0:taille:2] + Vz[1:taille:2]*1j;
  fileobj.close()
  A = pylab.zeros([n, n]) + 0j*pylab.zeros([n, n])
  ind  = 0
  for i in range(n):
    taille = n-i
    A[i,i:n] = V[ind:(ind+taille)];
    A[i:n,i] = V[ind:(ind+taille)];
    ind = ind + taille
  
  return A

def loadMat(nom_fichier):
  """ Reads a sparse matrix in coordinate format 
     Usage : A = loadMat(nom_fichier) """
  A = pylab.loadtxt(nom_fichier)
  B = scipy.sparse.coo_matrix((A[:,2],(pylab.array(A[:,0]-1, dtype='int64'),pylab.array(A[:,1]-1, dtype='int64'))))
  return scipy.sparse.csr_matrix(B)
  
def loadComplexMat(nom_fichier):
  """ Reads a complex sparse matrix in coordinate format 
     Usage : A = loadComplexMat(nom_fichier) """
  fid = open(nom_fichier,'r')
  V = fid.readlines()
  fid.close()
  n = len(V)
  val = pylab.zeros(n) + 1j*pylab.zeros(n)
  row = list(range(n))
  col = list(range(n))
  for i in range(n):
    mots = V[i].split()
    row[i] = int(mots[0])-1
    col[i] = int(mots[1])-1
    pos_comma = mots[2].find(',')
    real_part = float(mots[2][1:pos_comma])
    imag_part = float(mots[2][pos_comma+1:len(mots[2])-1])
    val[i] = real_part + 1j*imag_part
    
  A = scipy.sparse.coo_matrix((val,(row,col)))
  return scipy.sparse.csr_matrix(A)

def loadArray3D_real(nom_fichier, double_prec = True):
  """ Reads a 3-D array (Seldon structure) """
  fileobj = open(nom_fichier, mode='rb')
  m = pylab.fromfile(fileobj, 'i', 1)[0]
  n = pylab.fromfile(fileobj, 'i', 1)[0]
  k = pylab.fromfile(fileobj, 'i', 1)[0]
  taille = m*n*k
  if (double_prec):
    V = pylab.fromfile(fileobj, 'd', taille)
  else:
    V = pylab.fromfile(fileobj, 'f', taille).astype('float64')
  
  fileobj.close()
  A = pylab.reshape(V,(m, n, k))
  return A

def loadComplexVec(nom_fichier):
  """ Reads a complex vector in ascii format
     Usage : A = loadComplexVec(nom_fichier) """
  fid = open(nom_fichier,'r')
  V = fid.readlines()
  fid.close()
  n = len(V)
  val = pylab.zeros(n) + 1j*pylab.zeros(n)
  for i in range(n):
    ligne = V[i].strip();
    pos_comma = ligne.find(',')
    real_part = float(ligne[1:pos_comma])
    imag_part = float(ligne[pos_comma+1:len(ligne)-1])
    val[i] = real_part + 1j*imag_part
    
  return val

def trajectoire_racine(z, R, Ri, sauve_film = None):
  """ displays the evolution in time of roots
   R and Ri are real part and imaginary part of the roots
   for different times, they are matrices, each root being stored in a column
   z are the time values

  Usage : film1D(z, R, Ri)
          film1D(z, R, Ri, file_name)

  Put a file_name if you want to save the frames in png files """
  pylab.figure();
  pylab.axis('image')
  pylab.axis(pylab.array([-1.5, 1.5, -2, 2.0]))
  lineR = [0]*R.shape[1]
  markerR = [0]*R.shape[1]
  colorR = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
  colorR += colorR
  for j in range(R.shape[1]):
    lineR[j], = pylab.plot(R[0:2, j], Ri[0:2, j], colorR[j])
    markerR[j], = pylab.plot(R[1, j], Ri[1, j], colorR[j] + 'o')
    
  for i in range(3, R.shape[0]):
    for j in range(R.shape[1]):
      lineR[j].set_data(R[0:i, j], Ri[0:i, j])
      markerR[j].set_data(R[i-1:i, j], Ri[i-1:i, j])
    
    pylab.title("z = " + str(z[i]))
    pylab.draw()
    if (sauve_film != None):
      pylab.savefig(sauve_film+EntierToString(i)+".png")
  

def film1D(V, cmin = 'auto', cmax = 'auto', dt = 0, sauve_film = None):
  """ displays the evolution in time of a 1-D solution
   V is assumed to be a matrix, each column of the matrix contains
   the solution at a different time step

  Usage : film1D(V)
          film1D(V, dt)

  Put a dt different from 0, if the animation is too fast, so that the system
  will wait dt seconds between each display """
  N = pylab.size(V, 0);
  x = pylab.linspace(-1, 1, N);
  Vmin = float(min(pylab.reshape(V, (pylab.size(V), 1))));
  Vmax = float(max(pylab.reshape(V, (pylab.size(V), 1))));
  deltaV = Vmax - Vmin;
  if (cmin == 'auto'):
    Vmin = Vmin - 0.1*deltaV;
  else:
    Vmin = cmin
    
  if (cmax == 'auto'):
    Vmax = Vmax + 0.1*deltaV;
  else:
    Vmax = cmax
  
  line, = pylab.plot(x, V[:,0]);
  pylab.axis(pylab.array([-1, 1, Vmin, Vmax]));
  os.system('sleep ' +str(dt))
  for i in range(1, pylab.size(V, 1)):  
    line.set_ydata(V[:,i])
    pylab.draw()
    if (sauve_film != None):
      pylab.savefig(sauve_film+EntierToString(i)+".png")

    os.system('sleep ' +str(dt))

def film1D_file(base, ext, debut, fin, comp, vmin, vmax,
                fichier_teta = None, coef_ampli = 1.0, dt = 0, rapport_teta = 1,
                sauve_film = None):
  """ displays the evolution in time of a 1-D solution
   V is assumed to be stored in different files 
   of the form base0ext.dat, base1ext.dat, etc

  Usage : film1D_file("Val", "shank.txt", 0, 100, 0)
          film1D_file(base, ext, debut, fin, comp, vmin, vmax, fichier_teta, coef_ampli, dt, rapport_teta, sauve_film)
          
  base : base of the files containing the displacements of the shank
  ext : extension of the files containing the displacements of the shank
        The different files will be base0ext, base1ext, etc
        For instance, if base is equal to Val, ext to shank.txt,
        the files are Val0shank.txt, Val1shank.txt, etc
  debut : first snapshot number to display (eg. if equal to 23, we start from Val23shank.txt)
  fin : last snapshot number to display
  comp : component to display (0 for w, 1 for phi), i.e. column of Val?shank.txt considered
  vmin,vmax : scales displayed in y-axis
  fichier_teta : file where teta is stored (usually Sismo_Shank.txt)
                 if not provided, we consider that the shank is not rotating
  coef_ampli : amplification coefficient to apply to the displacement in order to exagerate the vibration
  dt : waiting time before displaying the next snapshot
  rapport_teta : if the value teta is written each dt_teta time, and the displacement each dt time
                 rapport_teta is equal to dt_teta/dt
  sauve_film : if provided, the snapshots are written on the disk, with a file name beginning with sauve_film
  
  Put a dt different from 0, if the animation is too fast, so that the system
  will wait dt seconds between each display """
  # si on a un teta en plus
  if (fichier_teta != None):
    teta = pylab.loadtxt(fichier_teta);
  
  # the last snapshot is loaded
  inst = fin-1
  if (comp == 'binary'):
    [X, Y, Z, coor, V] = loadND(base + EntierToString(inst) + ext);
  else:
    V = pylab.loadtxt(base + str(inst) + ext);
  
  N = V.shape[0]
  L = 0.13
  Hp = 0.01
  sL = pylab.linspace(0, L, N);
  sH = pylab.linspace(0, Hp, N);
  
  # 
  if (comp == 'binary'):
    [X, Y, Z, coor, V] = loadND(base + EntierToString(debut) + ext);
    Vc = pylab.real(V)
  else:
    V = pylab.loadtxt(base + str(debut) + ext);
    Vc = V[:,comp]
  
  if (fichier_teta != None):
    x = sL*cos(teta[0, 3]) - coef_ampli*(Vc*sin(teta[0,3]));
    y = sL*sin(teta[0, 3]) + coef_ampli*(Vc*cos(teta[0,3]));
    line, = pylab.plot(x, y)
    xH = L*cos(teta[0,3])-sH*sin(teta[0,3])
    yH = L*sin(teta[0,3])+sH*cos(teta[0,3])
    lineH, = pylab.plot(xH, yH)

    pylab.axis('equal')
    pylab.axis(pylab.array([-1.32*L, 1.32*L, -1.32*L, 1.32*L]));
  else:
    line, = pylab.plot(X, Vc)
    pylab.axis(pylab.array([0, 10.0, vmin, vmax]));  

  os.system('sleep ' +str(dt))
  for i in range(debut+1, fin):  
    if (comp == 'binary'):
      [X, Y, Z, coor, V] = loadND(base + EntierToString(i) + ext);
      Vc = pylab.real(V)
    else:
      V = pylab.loadtxt(base + str(i) + ext);
      Vc = V[:, comp]
    
    if (fichier_teta != None):
      x = sL*cos(teta[rapport_teta*i, 3]) - coef_ampli*(Vc*sin(teta[i*rapport_teta,3]));
      y = sL*sin(teta[rapport_teta*i, 3]) + coef_ampli*(Vc*cos(teta[i*rapport_teta,3]));
      line.set_xdata(x)
      line.set_ydata(y)
      xH = x[len(x)-1]-sH*sin(teta[rapport_teta*i,3])
      yH = y[len(y)-1]+sH*cos(teta[rapport_teta*i,3])
      lineH.set_xdata(xH)
      lineH.set_ydata(yH)      
    else:
      line.set_ydata(Vc)
    
    pylab.draw()
    t_courant = '{:04.3f}'.format(5e-4*i)
    pylab.title('t = ' + t_courant + '  x'+str(coef_ampli))
    if (sauve_film != None):
      pylab.savefig(sauve_film+EntierToString(i)+".png")
      
    if (dt != 0):
      os.system('sleep ' +str(dt))
  
  
def EntierToString(i):
  """ converts an integer into a string of four characters
  Usage : S = EntierToString(i)
   """
  if (i < 10):
    return '000'+str(i)
  elif (i < 100):
    return '00'+str(i)
  elif (i < 1000):
    return '0'+str(i)
  else:
    return str(i)
  
def affiche_mode(root, ext, n1, n2, inc=1, cplx=0):
  """ Displays 2-D solutions, which are stored in files 
     of the type toto0000Ext.dat, toto0001Ext.dat, toto0002Ext.dat 

  Usage : affiche_mode('toto', 'Ext', 0, 500) 
          affiche_mode(base, ext, n1, n2, inc)

     inc defines the increment, so that you can display only even numbers for example. """
  mode_to_keep = [];
  for i in range(n1, n2, inc):
    nom = root + EntierToString(i) + ext + ".dat"
    print(nom)
    [X, Y, Z, coor, V] = loadND(nom)
    if (len(V) == 3):
      plot3d_plane(X, Y, Z, V, coor)
    else:
      if (Y.min() == Y.max()):
        if (cplx == 0):
          plot2dinst(X, Z, V)
        else:
          plot2dinst(X, Z, pylab.real(V))
          plot2dinst(X, Z, pylab.imag(V))
      else:
        if (cplx == 0):
          plot2dinst(X, Y, V)
        else:
          plot2dinst(X, Y, pylab.real(V))
          plot2dinst(X, Y, pylab.imag(V))
    
    line = input("Press o to keep the mode")
    if (line == "o"):
      mode_to_keep.append(i)
      
  return mode_to_keep

def calcule_vitesse(t, U):
  """ Computes the velocity from the displacements
  Usage : V = calcule_vitesse(t, U)

  t : times (regular subdivision as if produced by arange)
  U : displacements 
  V : velocity, approximation of dU/dt
  
  The velocities are computed with a fourth-order approximation at the center,
  and a second-order approximation at the two extremities
  """
  dt = (max(t) - min(t)) / len(t)
  V = U.copy()
  nby = len(U)
  V[2:nby-2] = (-U[4:nby] + 8*U[3:nby-1] - 8*U[1:nby-3] + U[0:nby-4])/(12.0*dt)
  V[1] = (U[2] - U[0])/(2.0*dt);
  V[nby-2] = (U[nby-1] - U[nby-3])/(2.0*dt);
  V[0] = (-3.0*U[0] + 4.0*U[1] - U[2])/(2.0*dt);
  V[nby-1] = (3.0*U[nby-1] - 4.0*U[nby-2] + U[nby-3])/(2.0*dt);
  return V

def calcule_acc(t, U):
  """ Computes the acceleration from the displacements
  Usage : A = calcule_acc(t, U)

  t : times (regular subdivision as if produced by arange)
  U : displacements 
  A : acceleration, approximation of d^2 U/dt^2
  
  calcule_vitesse is called twice to obtain the acceleration """
  V = calcule_vitesse(t, U)
  A = calcule_vitesse(t, V)
  return A

def wavread(file_name):
  """ reads a .wav file 
      
      Usage : y = wavread('son.wav') 
      
      file_name : nom du fichier .wav a lire 
      y : donnees lues (flottants)  """
    
  fid = wave.open(file_name, 'r')
  
  # nombre de frames
  nframes = fid.getnframes()
  
  # nombre de bytes par echantillon
  nbytes = fid.getsampwidth()
  
  # frequence d'echantillonage
  nfreq = fid.getframerate()
  
  # donnees
  data = fid.readframes(nframes)

  # nombre de channels
  nchannel = fid.getnchannels()
  
  fid.close()
  
  print("Frequence d'echantillonage", nfreq)

  if (nbytes == 2):
    v = struct.unpack('h'*nframes*nchannel, data)
    v = pylab.array(v)/32767.0
    V = pylab.zeros([nframes, nchannel])
    for n in range(nchannel):
      V[:,n] = v[n:nframes*nchannel:nchannel]
    
    return V
  else:
    print("Not implemented")

def wavwrite(y, file_name, nfreq, vmax = 0):
  """ Writes a .wav file
  
  Usage : wavwrite(data, "son.wav", 48000) 
  
  y : signal a ecrire
  file_name : nom du fichier .wav
  nfreq : frequence d'echantillonage
  vmax (optionnel) : amplitude maximum du signal """
  
  fid = wave.open(file_name, 'w')
  
  fid.setsampwidth(2)
  fid.setframerate(nfreq)
  
  if (len(y.shape) == 1):
    nframe = len(y)
    nchannel = 1
    fid.setnframes(nframe)
    if (vmax == 0):
      vmax = 1.01*max(vmax, max(abs(y)))
    
    fid.setnchannels(1)
    v = y / vmax * 32767
    num = [0]*len(v)
    for i in range(len(v)):
      num[i] = int(v[i])
  else:
    nframe = y.shape[0]
    fid.setnframes(nframe)
    nchannel = y.shape[1]
    if (vmax == 0):
      vmax = 1.01*max(vmax, abs(y).max())
    
    fid.setnchannels(nchannel)
    num = [0]*nchannel*nframe
    for n in range(nchannel):
      for i in range(y.shape[0]):
        v = y[i,n] / vmax * 32767
        num[n+i*nchannel] = int(v)
  
  data = struct.pack('h'*nframe*nchannel, *num)
  fid.writeframes(data)
  fid.close()

def ReadSismo(nom_fichier, MAXI = 1e7, no_time = False):
  """ Reading of a "Montjoie" seismogramm in binary format
  usage : sis = ReadSismo(nom_fichier)
          sis = ReadSismo(nom_fichier, MAXI)
 
  nom_fichier : name of the output file
  MAXI : maximum number of entries to read in the output file
  sis : seismogramm, first column of sis contains the time, other columns values """
  
  # opening the file
  fid = open(nom_fichier, mode='rb')
  
  # type of data (0 : float, 1 : double)
  type_data = pylab.fromfile(fid, 'i', 1)[0]
  prec = 'd'
  if (type_data == 0):
    prec = 'f';
    
  compt = 0
  while (compt < MAXI):
    try:
      n = pylab.fromfile(fid, 'i', 1)[0]
      Vloc = pylab.fromfile(fid, prec, n)
      if (compt == 0):
        L = pylab.zeros([MAXI, len(Vloc)])
      
      L[compt, :] = Vloc;   
      compt = compt+1;
    except:
      # fin du fichier
      break
    
  if (no_time):
    L = L[0:compt, :];
  else:
    L = L[0:compt, :];
    V = L;
    b = GetUniqueIndex(V[:,0])
    L = V[b,:]
  
  return L

def ReadFarField(name_file, MAXI):
  """ Reads a far field produced by Montjoie
  usage : Points, V = ReadFarField(name_file, nb_max) 
  
  name_file : name where the far field is stored
  nb_max : maximum number of values to read """
  
  fid = open(name_file,mode='rb')

  # type of data (0 : float, 1 : double, 2 : complex float, 3 : complex double)
  type_data = pylab.fromfile(fid, 'i', 1)[0]

  prec = 'd'; prec_size = 8;
  if ((type_data == 0) or (type_data == 2)):
    prec = 'f'
    prec_size = 4;
    
  n = pylab.fromfile(fid, 'i', 1)[0]
  Points = pylab.fromfile(fid, prec, n)
  Points = pylab.reshape(Points, [3,n/3])

  L = pylab.zeros([MAXI, n/3+1]);
  
  compt=0;
  while(compt<MAXI):
    try:
      n = pylab.fromfile(fid, 'i', 1)[0]
      Vloc = pylab.fromfile(fid, prec, n)
      L[compt, :] = Vloc;   
      compt=compt+1;
    except:
      # fin du fichier
      break
    
  L = L[0:compt, :];
  V = L;
  b = GetUniqueIndex(V[:,0])
  L = V[b,:]
  return Points, L

def ReadNodalPoints(fichier):
  """ reads nodal points from a given file 
  usage : pts, nodes = ReadNodalPoints(file_name) """
  fid = open(fichier, 'r')
  lignes = fid.readlines()
  mots = lignes[0].split(' ')
  nb_elt = int(mots[0])
  nb_pts = int(mots[1])
  pts = pylab.loadtxt(lignes[1:nb_pts+1])
  nodes = [0]*nb_elt
  for i in range(nb_elt):
    nodes[i] = pylab.array(lignes[nb_pts+1+2*i+1].split(' '), dtype='int32')
  
  return pts, nodes

def ReadMeshData(fichier):
  """ reads an index from a given file
  usage : nu = ReadMeshData(fichier) """
  fid = open(fichier, "r");
  type_data = pylab.fromfile(fid, 'i', 1)[0]
  print("type_data = ", type_data)
  nb_elt = pylab.fromfile(fid, 'i', 1)[0]
  offset = pylab.fromfile(fid, 'i', nb_elt+1)
  nu = [0]*nb_elt
  if (type_data == 0):
    for i in range(nb_elt):
      nu[i] = pylab.fromfile(fid, 'd', offset[i+1]-offset[i])
  else:
    for i in range(nb_elt):
      taille = 2*(offset[i+1]-offset[i])
      data = pylab.fromfile(fid, 'd', taille)
      nu[i] = data[0:taille:2] + 1j*data[1:taille:2]
  
  return nu

def WriteMeshData(nu, fichier):
  """ writes an index into a file 
  usage : WriteMeshData(nu, fichier) """
  fid = open(fichier, "w")
  nb_elt = len(nu)
  if (nu[0].dtype == 'float64'):
    type_data = pylab.array([0, nb_elt], dtype='int32')
  elif (nu[0].dtype == 'complex128'):
    type_data = pylab.array([1, nb_elt], dtype='int32')
  
  type_data.tofile(fid)
  offset = pylab.array([0]*(nb_elt+1)).astype('int32')
  for i in range(nb_elt):
    offset[i+1] = offset[i] + len(nu[i])
    
  offset.tofile(fid)
  for i in range(nb_elt):
    nu[i].tofile(fid)
  
  fid.close()

def spectrocoupe(y, nfft, fs, npad, nover, seuil_min = -60):
  """ Displays a spectrogram of a data  in decibels
  Usage : spectrocoupe(y, nfft, fs, npad, nover) 
          spectrocoupe(y, nfft, fs, npad, nover, min_threshold)
          
          y : data to analyse
          nfft : number of points to use for computation of fft
          fs : sample frequency of the data
          npad : extra-zeros to add in order to perform the fft
          nover : number of points for the overlap between two ffts
          min_threshold : by default equal to -60 dB """
  
  fig = pylab.figure()
  [B, f, t, im] = pylab.specgram(y, nfft, fs, pad_to=npad, noverlap=nover, hold=False)
  #pylab.close()
  bmin = B.max()
  fmin = 0
  fmax = 0.5*fs
  tmin = t[0]
  tmax = t[len(t)-1]  
  pylab.imshow(20.0*pylab.log10(abs(B)/bmin), vmin = seuil_min, vmax = 0, interpolation='bilinear', \
                 origin='lower', extent=(tmin, tmax, fmin, fmax), aspect='auto',cmap=pylab.cm.jet)
  
  pylab.ylim(0, fmax/2.5)
  pylab.colorbar()
  pylab.xlabel('Temps (s)');
  pylab.ylabel('Frequence (Hz)');

def GetFourier(V, dt):
  """ Computes Fourier transform of V
  Usage : omega, Vhat = GetFourier(V, dt) 
  
          dt : time step between each value of V
          omega : pulsations for which Fourier transform is computed
          Vhat : Fourier transform for each point in omega  """
  
  coef_fft = dt / pylab.sqrt(2.0*pi)
  Vchap = coef_fft*pylab.fft(V)
  omega = pylab.linspace(-pi/dt, pi/dt, len(V)+1)
  omega = omega[0:len(V)]
  return omega, pylab.fftshift(Vchap)

def trace_fft(signal, F):
  """ Displays Fourier transform of a signal in decibels
   Usage : trace_fft(y, f)
   
   y : data to analyse
   f : sample frequency of the data  """
  
  fft_sig = pylab.fft(signal);
  max_fft = max(abs(fft_sig));
  nfreq = len(signal)
  freq = pylab.linspace(-0.5*F, 0.5*F, len(signal)+1)
  pylab.plot(freq[0:nfreq], pylab.fftshift(20.0*pylab.log10(abs(fft_sig)/max_fft)));
  pylab.xlim(0,0.4*F);
  pylab.xlabel('Frequency');
  pylab.ylabel('Fourier Transform : log(abs)')

def deriveX(X, Y, V):
  """ Computes the derivative of a 2-D field with respect to x-coordinate
  Usage : dVx = deriveX(X, Y, V)
  
          X : regularly spaced x-coordinate (as produced by linspace or meshgrid)
          Y : regularly spaced y-coordinate (as produced by linspace or meshgrid)
          V : 2-D field associated with points (X, Y) 
          dVx : derivative of V with respect to X
          The derivative is computed with fourth-order finite-differences inside the rectangle,
          and with second-order finite differences for points on the boundary"""       
  xmin = X.min()
  xmax = X.max()
  nbx = V.shape[0]
  nby = V.shape[1]
  dx = (xmax-xmin)/(nby-1);
  U = V.copy()
  U[:, 2:nby-2] = (-V[:, 4:nby] + 8*V[:,3:nby-1] - 8*V[:, 1:nby-3] + V[:, 0:nby-4])/(12.0*dx);
  
  U[:, 1] = (V[:, 2] - V[:, 0])/(2.0*dx);     
  U[:, nby-2] = (V[:, nby-1] - V[:, nby-3])/(2.0*dx);
  U[:, 0] = (-3.0*V[:, 0] + 4.0*V[:, 1] - V[:, 2])/(2.0*dx);
  U[:, nby-1] = (3.0*V[:, nby-1] - 4.0*V[:, nby-2] + V[:, nby-3])/(2.0*dx);
  return U

def deriveY(X, Y, V):
  """ Computes the derivative of a 2-D field with respect to y-coordinate
  Usage : dVy = deriveY(X, Y, V)
  
          X : regularly spaced x-coordinate (as produced by linspace or meshgrid)
          Y : regularly spaced y-coordinate (as produced by linspace or meshgrid)
          V : 2-D field associated with points (X, Y) 
          dVy : derivative of V with respect to Y
          The derivative is computed with fourth-order finite-differences inside the rectangle,
          and with second-order finite differences for points on the boundary"""       
  ymin = Y.min()
  ymax = Y.max()
  nbx = V.shape[1]
  nby = V.shape[0]
  dy = (ymax-ymin)/(nby-1);
  U = V.copy()
  U[2:nby-2, :] = (-V[4:nby, :] + 8*V[3:nby-1, :] - 8*V[1:nby-3, :] + V[0:nby-4, :])/(12.0*dy);
  
  U[1, :] = (V[2, :] - V[0, :])/(2.0*dy);     
  U[nby-2, :] = (V[nby-1, :] - V[nby-3, :])/(2.0*dy);
  U[0, :] = (-3.0*V[0, :] + 4.0*V[1, :] - V[2, :])/(2.0*dy);
  U[nby-1, :] = (3.0*V[nby-1, :] - 4.0*V[nby-2, :] + V[nby-3, :])/(2.0*dy);
  return U

def log_triangle(pos, order, taille, couleur, c0_x = 0.35, c0_y=0.1, c1_x = -0.3, c1_y = -0.6):
  """ Displays a small triangle (indicating the order in a convergence plot)
  Usage : log_triangle(pos, order, size, color)
          log_triangle(pos, order, size, color, c0_x, c0_y, c1_x, c1_y)

          pos : position where the triangle should be plotted (eg. [1e-2, 1e-4])
          order : order of convergence to display (eg. 4)
          size : size of the triangle (eg. 0.12)
          color : color of the triangle (eg. 'm' for magenta)
          c0_x, c0_y : position of the text containing '1' (default : 0.35 and 0.1)
          c1_x, c1_y : position of the text containing '4' (default : -0.3 and -0.6) """
  x0 = pos[0]
  y0 = pos[1]
  dx = taille
  dy = taille*order
  pylab.plot([x0, x0*10**dx, x0, x0], [y0, y0, y0*10**(-dy), y0], couleur)
  pylab.text(x0*10**(c0_x*dx), y0*10**(c0_y*dy), '1')
  pylab.text(x0*10**(c1_x*dx), y0*10**(c1_y*dy), str(abs(order)))

def ReadPointsMesh(nom):
  """ Reads points contained in a .mesh file
  Usage : x, y = ReadPointsMesh(file_name)
  
          file_name : name of the .mesh file
          x : x-coordinates of the points stored in the file
          y : y-coordinates of the points stored in the file """
  fid = open(nom, "r")
  test_loop = True
  vertices_found = False
  while (test_loop):
    ligne = fid.readline()
    if (len(ligne) == 0):
      test_loop = False
    else:
      if (ligne[len(ligne)-1] != '\n'):
        test_loop = False
      
    if (len(ligne) > 2):
      if (ligne[0:len(ligne)-1].strip() == 'Vertices'):
        test_loop = False
        vertices_found = True
        
  if (vertices_found):
    nb_vert = int(fid.readline())
    #print("nb vertices" , nb_vert)
    x = pylab.zeros(nb_vert)
    y = pylab.zeros(nb_vert)
    for i in range(nb_vert):
      ligne = fid.readline()
      mots = ligne.split()
      x[i] = float(mots[0])
      y[i] = float(mots[1])
      
    return x, y

def ReadBb(nom):
  """ Reads a .bb file
  Usage : u = ReadBb(file_name)
  
          file_name : name of the file .bb
          u : solution stored in the file """
  fid = open(nom, 'r')
  s = fid.readlines()
  if (len(s) > 1):
    t = pylab.array(s[1:])
    return t.astype(float)

def ReadComplexBb(nom):
  """ Reads a complex solution stored in .bb files
  Usage : u = ReadBb(file_name)
  
          file_name : solution is stored in file_name_real.bb and file_name_imag.bb
          u : solution stored in the two files ending with _real.bb and _imag.bb """
  return ReadBb(nom+'_real.bb') + 1j*ReadBb(nom+'_imag.bb')

def loadComplexMatPar(fic, n, assemble = True):
  m = 0;
  for i in range(n):
    A = loadComplexMat(fic + '_P' + str(i) + '.dat');
    m = max(m, pylab.size(A, 0));
    m = max(m, pylab.size(A, 1));
    
  
  A = pylab.zeros([m, m]) + 1j*pylab.zeros([m, m])
  for i in range(n):
    B = loadComplexMat(fic + '_P' + str(i) + '.dat').todense();
    if (assemble):
      A[0:pylab.size(B,0), :][:, 0:pylab.size(B,1)] += B;
    else:
      A[0:pylab.size(B,0), :][:, 0:pylab.size(B,1)] = B;

  return A

def loadMatPar(fic, n, assemble = True):
  m = 0;
  for i in range(n):
    A = loadMat(fic + '_P' + str(i) + '.dat');
    m = max(m, pylab.size(A, 0));
    m = max(m, pylab.size(A, 1));
    
  
  A = pylab.zeros([m, m]) + 1j*pylab.zeros([m, m])
  for i in range(n):
    B = pylab.array(loadMat(fic + '_P' + str(i) + '.dat').todense());
    if (assemble):
      A[0:pylab.size(B,0), 0:pylab.size(B,1)] += B;
    else:
      A[0:pylab.size(B,0), 0:pylab.size(B,1)] = B;

  return A

def write_txtFile(fic,Sol):
  tab = Sol.tolist();
  lines_to_write=[];
  for k in range(len(tab)):
    real_part = tab[k].real;
    imag_part = tab[k].imag;
    lines_to_write.append('['+ str(real_part) + ',' + str(imag_part) +']\n');
  file = open('../../../Documents/THESE/out/' + fic + '.txt','w');
  file.writelines(lines_to_write);
  file.close();

def write_binaryFile(fic,Sol,OutputPath):
  tab = Sol.tolist();
  
  M = np.zeros(2*len(tab));
  for p in range(len(tab)-1):
    M[2*p] = np.real(tab[p]);
    M[2*p+1] = np.imag(tab[p]);

  file = open(OutputPath + fic, 'w');
  M.astype('float64').tofile(file)
  file.close();
  
  return M

def write_binaryDirectory(InputPath,OutputPath):
  compt = 0;
  list_elements =  os.listdir(InputPath);
  list_elements.sort();
  for element in list_elements:
    if element.endswith('.dat'):
      rest = compt % 3;
      ntest = int(np.floor(compt/3 + 1));
      [X,Y,Z,coor,u] = loadND(InputPath + element);
      if (rest == 0):
        fic = 'tab' + str(ntest) + '_MONTJOIE_x.dat';
        write_binaryFile(fic,u,OutputPath);
      elif (rest == 1):
        fic = 'tab' + str(ntest) + '_MONTJOIE_y.dat';
        write_binaryFile(fic,u,OutputPath);
      else:
        fic = 'tab' + str(ntest) + '_MONTJOIE_z.dat';
        write_binaryFile(fic,u,OutputPath);
      compt = compt+1;  
  return compt
  
def write_binaryRCS(fic,Sol,OutputPath):
  tab = Sol.tolist();
  print(tab[0]);
  
  M = np.zeros(2*len(tab));
  for p in range(len(tab)):
    M[2*p] = tab[p][0];
    M[2*p+1] = tab[p][1];
  print(M[0]);
  print(M[1]);
  file = open(OutputPath + fic ,'w');
  M.astype('float32').tofile(file);
  file.close();
  return M

def plot_eigenval(N):
  list_mark = ['+','.','o','v','<','>','s','p','h','D']
  list_col = ['b','g','r','c','m','y','k']
  list_sign = [0]*len(list_mark)*len(list_col)
  for i in range(len(list_mark)):
    for j in range(len(list_col)):
      list_sign[j*len(list_mark) + i] = list_mark[i] + list_col[j]

  permut = np.random.permutation(70)
  for i in range(N):
    if (i != 21):
      L = load1D('../EigenvalSphereS' + str(i) + '__sorted.dat')
      i0 = permut[i%70]
      if ((list_sign[i0][0] == "+") or (list_sign[i0][0] == ".")):
        pylab.plot(-np.imag(L), -np.real(L), list_sign[i0])
      else:
        pylab.plot(-np.imag(L), -np.real(L), list_sign[i0], fillstyle='none')

def plot_circle(nom_fichier):
  fig = pylab.figure()
  ax = fig.gca();
  shift = pylab.loadtxt(nom_fichier)
  list_col = ['b','g','r','c','m','y','k']
  N = len(list_col)
  Nshift = shift.shape[0]
  xmin = 1e300; xmax = -1e300; ymin = 1e300; ymax = -1e300;
  for i in range(Nshift):
    R = shift[i, 2]
    omegaR = shift[i, 0]
    omegaI = shift[i, 1]
    xmin = min(xmin, omegaR-R)
    xmax = max(xmax, omegaR+R)
    ymin = min(ymin, omegaI-R)
    ymax = max(ymax, omegaI+R)
    c = pylab.Circle((omegaR, omegaI), R, color=list_col[i%N], fill=False)
    ax.add_artist(c)

  pylab.axis('equal')
  pylab.xlim(xmin, xmax)
  pylab.ylim(ymin, ymax)
  
def load_optical_index(nom):
  fid = open(nom)
  lignes = fid.readlines();
  s = lignes[2].split("*");
  N = int(s[1])
  x = pylab.zeros(N); y = pylab.zeros(N); z = pylab.zeros(N)
  for k in range(N):
    s = lignes[k+3].split("*")
    x[k] = float(s[2])
    y[k] = float(s[3])
    z[k] = float(s[4])

  return x, y, z

def count_eigenval(L, om1_r, om2_r, om1_i, om2_i):
  n = 0
  for i in range(len(L)):
    if ((pylab.real(L[i]) < om2_r) and (pylab.real(L[i]) > om1_r) and (pylab.imag(L[i]) < om2_i) and (pylab.imag(L[i]) > om1_i)):
      n += 1

  return n

def getResidualFile(name):
  fid = open(name)
  lignes = fid.readlines()
  n = []
  res = []
  num_prev = -1
  for i in range(len(lignes)):
    if (lignes[i].startswith('Residu at iteration number')):
      mots = lignes[i][27:].split()
      num = int(mots[0])
      if (num != num_prev):
        n.append(num)
        res.append(float(mots[1]))

      num_prev = num

  return np.array(n), np.array(res)
