#!/usr/bin/env python 
""" Module for reading Montjoie output files, and displaying them,
   Python equivalent of Matlab scripts"""

import numpy as np
import scipy
import scipy.sparse
import os, sys, pickle
#import pylab
import struct, wave
#from pylab import cos, sin, pi
from scipy.interpolate import RectBivariateSpline
from scipy.interpolate import interp1d

def Refine1D(xmin, xmax, U, raff):
  nbx = len(U)
  x = np.linspace(xmin, xmax, nbx);
  xr = np.linspace(xmin, xmax, raff*(nbx-1)+1);
  Up = interp1d(x, U, 'cubic')
  return Up(xr)

def Refine2D(xmin, xmax, ymin, ymax, U, raff):
  nbx = U.shape[0]
  nby = U.shape[1]
  x = np.linspace(xmin, xmax, nbx);
  y = np.linspace(ymin, ymax, nby);
  xr = np.linspace(xmin, xmax, raff*(nbx-1)+1);
  yr = np.linspace(ymin, ymax, raff*(nby-1)+1);
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
  uchap = np.fft.fft(u[0:N-1])
  uchap_pad = np.zeros(s-1) + 1j*np.zeros(s-1)
  # we fill uchap_pad with additional zeros (for high frequency)
  uchap_pad[0:(N-1)//2] = uchap[0:(N-1)//2]
  uchap_pad[s-1-(N-1)//2:] = uchap[(N-1)//2:]
  ur = np.zeros(s) + 1j*np.zeros(s)
  # we perform inverse fft to obtain ur on refined points
  ur[0:s-1] = (s-1) / float(N-1) * np.fft.ifft(uchap_pad)
  # periodic condition
  ur[s-1] = ur[0]
  # we return only real part (u is assumed to be a real vector)
  return np.real(ur)

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
  nb_grids = np.fromfile(fileobj, 'i', 1)[0]
  
  # the dimension
  dim = np.fromfile(fileobj, 'i', 1)[0]
  
  # type of data (0 : float, 1 : double, 2 : complex float, 3 : complex double)
  type_data = np.fromfile(fileobj, 'i', 1)[0]
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
  type_grid = np.fromfile(fileobj, 'i', 1)
  
  if (dim == 1):
    xmin = np.fromfile(fileobj, prec, 1)[0]
    xmax = np.fromfile(fileobj, prec, 1)[0]
    nbx = np.fromfile(fileobj, 'i', 1)[0]
    n = np.fromfile(fileobj, 'i', 1)[0]
    taille = int(cplx*nbx);
    data = np.fromfile(fileobj, prec, taille)
    if (cplx == 2):
      M = data[0:taille:2] + 1j*data[1:taille:2];
    else:
      M = data;
      
    fileobj.close()
    V = M
    if (raff != 1):
      V = Refine1D(xmin, xmax, V, raff)
      
    XI = np.linspace(xmin, xmax, (nbx-1)*raff + 1)
  elif (dim == 2):
    
    # depending the type of grid, we read the appropriate datas
    if (type_grid == 0):

      # SismoPlane
      offset = 0;
      xmin_ = np.zeros(nb_grids); xmax_ = np.zeros(nb_grids);
      ymin_ = np.zeros(nb_grids); ymax_ = np.zeros(nb_grids);
      nbx_ = np.zeros(nb_grids).astype('int'); nby_ = np.zeros(nb_grids).astype('int');
      for i in range(nb_grids):
        xmin_[i] = np.fromfile(fileobj, prec, 1)
        xmax_[i] = np.fromfile(fileobj, prec, 1)
        nbx_[i] = np.fromfile(fileobj, 'i', 1)[0]
        ymin_[i] = np.fromfile(fileobj, prec, 1)
        ymax_[i] = np.fromfile(fileobj, prec, 1)
        nby_[i] = np.fromfile(fileobj, 'i', 1)[0]     
        if (i < num_inst-1):
          offset = offset + nbx_[i]*nby_[i];
        
       
      xmin = xmin_[num_inst-1]; xmax = xmax_[num_inst-1];
      ymin = ymin_[num_inst-1]; ymax = ymax_[num_inst-1];
      nbx = nbx_[num_inst-1]; nby = nby_[num_inst-1];
      n = np.fromfile(fileobj, 'i', 1)[0]
      taille = int(cplx*nbx*nby);
      # we read values related to the required snapshot
      fileobj.seek(offset*prec_size, 1);
      data = np.fromfile(fileobj, prec, taille)
      if (cplx == 2):
        M = data[0:taille:2] + 1j*data[1:taille:2];
      else:
        M = data;
      
      fileobj.close()

      # here we obtain the transpose of the data written by the user
      # it is because 2-D displays naturally transpose the matrix (rows are displayed along y)
      U = np.reshape(M, (nby, nbx));
      if (raff != 1):
        U = Refine2D(ymin, ymax, xmin, xmax, U, raff)
      
      X = np.linspace(xmin, xmax, (nbx-1)*raff+1);
      Y = np.linspace(ymin, ymax, (nby-1)*raff+1);
      V = U;
      XI = X; YI = Y;
              
    elif (type_grid == 1):
      
      # SismoLine
      offset = 0;
      xmin_ = np.zeros(nb_grids); xmax_ = np.zeros(nb_grids);
      ymin_ = np.zeros(nb_grids); ymax_ = np.zeros(nb_grids);
      nbx_ = np.zeros(nb_grids).astype('int');
      for i in range(nb_grids):
        xmin_[i] = np.fromfile(fileobj, prec, 1)[0]
        xmax_[i] = np.fromfile(fileobj, prec, 1)[0]
        ymin_[i] = np.fromfile(fileobj, prec, 1)[0]
        ymax_[i] = np.fromfile(fileobj, prec, 1)[0]
        nbx_[i] = np.fromfile(fileobj, 'i', 1)[0]
        if (i < num_inst-1):
          offset = offset + nbx_[i];
        
      
      xmin = xmin_[num_inst-1]; xmax = xmax_[num_inst-1];
      ymin = ymin_[num_inst-1]; ymax = ymax_[num_inst-1];
      nbx = nbx_[num_inst-1];
      n = np.fromfile(fileobj, 'i', 1)[0]
      taille = int(cplx*nbx);
      # we read values related to the required snapshot
      fileobj.seek(offset*prec_size, 1);
      data = np.fromfile(fileobj, prec, taille)
      if (cplx == 2):
        M = data[0:taille:2] + 1j*data[1:taille:2];
      else:
        M = data;
      
      fileobj.close()
      
      U = M;
      V = U;
      if (raff != 1):
        V = Refine1D(0.0, 1.0, U, raff)

      X = np.linspace(xmin, xmax, (nbx-1)*raff+1);
      Y = np.linspace(ymin, ymax, (nbx-1)*raff+1);
      XI = X; YI = Y;
      
    elif (type_grid == 2):
      # SismoPoint
      offset = 0;
      xmin_ = np.zeros(nb_grids);
      ymin_ = np.zeros(nb_grids);
      for i in range(nb_grids):
        xmin_[i] = np.fromfile(fileobj, prec, 1)[0]
        ymin_[i] = np.fromfile(fileobj, prec, 1)[0]
        if (i < num_inst-1):
          offset = offset + 1
        
      xmin = xmin_[num_inst-1];
      ymin = ymin_[num_inst-1];
      n = np.fromfile(fileobj, 'i', 1)[0]
      taille = cplx;
      # we read values related to the required snapshot
      fileobj.seek(offset*prec_size, 1);
      data = np.fromfile(fileobj, prec, taille)
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
      xmin_ = np.zeros(nb_grids);
      ymin_ = np.zeros(nb_grids);
      radius_ = np.zeros(nb_grids);
      nbx_ = np.zeros(nb_grids).astype('int');
      for i in range(nb_grids):
        xmin_[i] = np.fromfile(fileobj, prec, 1)[0]
        ymin_[i] = np.fromfile(fileobj, prec, 1)[0]
        radius_[i] = np.fromfile(fileobj, prec, 1)[0]
        nbx_[i] = np.fromfile(fileobj, 'i', 1)[0]
        if (i < num_inst-1):
          offset = offset + nbx_[i];
        
      
      xmin = xmin_[num_inst-1];
      ymin = ymin_[num_inst-1];
      radius = radius_[num_inst-1];
      nbx = nbx_[num_inst-1];
      n = np.fromfile(fileobj, 'i', 1)[0]
      taille = int(cplx*nbx);
      # we read values related to the required snapshot
      fileobj.seek(offset*prec_size, 1);
      data = np.fromfile(fileobj, prec, taille)
      if (cplx == 2):
        M = data[0:taille:2] + 1j*data[1:taille:2];
      else:
        M = data;
      
      fileobj.close()
      
      teta = np.linspace(0.0, 2.0*np.pi, (nbx-1)*raff+1);
      X = radius*cos(teta) + xmin
      Y = radius*sin(teta) + ymin
      U = M
      V = U;
      if (raff != 1):
        V = Refine1D(0.0, 2.0*np.pi, U, raff)
      
      XI = X; YI = Y;
      
    elif (type_grid == 6):
      
      # SismoPointsFile
      offset = 0;
      nbx_ = np.zeros(nb_grids, dtype='int');
      for i in range(nb_grids):
        nbx_[i] = np.fromfile(fileobj, 'i', 1)[0]
        coord = np.fromfile(fileobj, prec, 2*int(nbx_[i]))
        if (i == num_inst-1):
          pts = coord
        
        if (i < num_inst-1):
          offset = offset + nbx_[i];          
        
      nbx = int(nbx_[num_inst-1]);
      n = np.fromfile(fileobj, 'i', 1)[0]
      taille = int(cplx*nbx);
      # we read values related to the required snapshot
      fileobj.seek(offset*prec_size, 1);
      data = np.fromfile(fileobj, prec, taille)
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
      xmin_ = np.zeros(nb_grids); xmax_ = np.zeros(nb_grids)
      ymin_ = np.zeros(nb_grids); ymax_ = np.zeros(nb_grids)
      zmin_ = np.zeros(nb_grids); zmax_ = np.zeros(nb_grids)
      nbx_ = range(nb_grids); nby_ = range(nb_grids); nbz_ = range(nb_grids);
      for i in range(nb_grids):
        xmin_[i] = np.fromfile(fileobj, prec, 1)[0]
        xmax_[i] = np.fromfile(fileobj, prec, 1)[0]
        nbx_[i] = np.fromfile(fileobj, 'i', 1)[0]
        ymin_[i] = np.fromfile(fileobj, prec, 1)[0]
        ymax_[i] = np.fromfile(fileobj, prec, 1)[0]
        nby_[i] = np.fromfile(fileobj, 'i', 1)[0]
        zmin_[i] = np.fromfile(fileobj, prec, 1)[0]
        zmax_[i] = np.fromfile(fileobj, prec, 1)[0]
        nbz_[i] = np.fromfile(fileobj, 'i', 1)[0]
        if (i < num_inst-1):
          offset = offset + nbx_[i]*nby_[i]*nbz_[i];
      
      n = np.fromfile(fileobj, 'i', 1)[0]
      xmin = xmin_[num_inst-1]; xmax = xmax_[num_inst-1]; nbx = nbx_[num_inst-1];
      ymin = ymin_[num_inst-1]; ymax = ymax_[num_inst-1]; nby = nby_[num_inst-1];
      zmin = zmin_[num_inst-1]; zmax = zmax_[num_inst-1]; nbz = nbz_[num_inst-1];
      taille = cplx*(nbx*nby*nbz)
      # we read values related to the required snapshot
      fileobj.seek(offset*prec_size, 1);
      data = np.fromfile(fileobj, prec, taille)
      if (cplx == 2):
        M = data[0:taille:2] + 1j*data[1:taille:2];
      else:
        M = data;
      
      fileobj.close()
      
      X = np.linspace(xmin, xmax, nbx);
      Y = np.linspace(ymin, ymax, nby);
      Z = np.linspace(zmin, zmax, nbz);
      if (raff != 1):
        print("Refinement not implemented in 3-D")
      
      # here we obtain a 3-D array with inverted x and z
      V = np.reshape(M, [nbz, nby, nbx]);
      XI = X; YI = Y; ZI = Z;
      
    elif (type_grid == 1):
      # SismoGrille
      offset = 0;
      xmin_ = np.zeros(nb_grids); xmax_ = np.zeros(nb_grids);
      ymin_ = np.zeros(nb_grids); ymax_ = np.zeros(nb_grids);
      zmin_ = np.zeros(nb_grids); zmax_ = np.zeros(nb_grids);
      nbx_ = np.zeros(nb_grids, dtype='int');
      nby_ = np.zeros(nb_grids, dtype='int');
      nbz_ = np.zeros(nb_grids, dtype='int');
      coorx_ = np.zeros(nb_grids);
      coory_ = np.zeros(nb_grids);
      coorz_ = np.zeros(nb_grids);
      for i in range(nb_grids):
        xmin_[i] = np.fromfile(fileobj, prec, 1)[0]
        xmax_[i] = np.fromfile(fileobj, prec, 1)[0]
        nbx_[i] = np.fromfile(fileobj, 'i', 1)[0]
        ymin_[i] = np.fromfile(fileobj, prec, 1)[0]
        ymax_[i] = np.fromfile(fileobj, prec, 1)[0]
        nby_[i] = np.fromfile(fileobj, 'i', 1)[0]
        zmin_[i] = np.fromfile(fileobj, prec, 1)[0]
        zmax_[i] = np.fromfile(fileobj, prec, 1)[0]
        nbz_[i] = np.fromfile(fileobj, 'i', 1)[0]
        coorx_[i] = np.fromfile(fileobj, prec, 1)[0]
        coory_[i] = np.fromfile(fileobj, prec, 1)[0]
        coorz_[i] = np.fromfile(fileobj, prec, 1)[0]
        if (i < num_inst-1):
          offset = offset + nbx_[i]*nby_[i] + nbx_[i]*nbz_[i] + nby_[i]*nbz_[i];
      
      xmin = xmin_[num_inst-1]; xmax = xmax_[num_inst-1];
      ymin = ymin_[num_inst-1]; ymax = ymax_[num_inst-1];
      zmin = zmin_[num_inst-1]; zmax = zmax_[num_inst-1];
      coor = np.zeros(3); coor[0] = coorx_[num_inst-1];
      coor[1] = coory_[num_inst-1]; coor[2] = coorz_[num_inst-1];
      nbx = nbx_[num_inst-1]; nby = nby_[num_inst-1]; nbz = nbz_[num_inst-1];
      n = np.fromfile(fileobj, 'i', 1)[0]
      taille = int(cplx*(nbx*nby + nbx*nbz + nby*nbz));
      # we read values related to the required snapshot
      fileobj.seek(offset*prec_size, 1);
      data = np.fromfile(fileobj, prec, taille)
      if (cplx == 2):
        M = data[0:taille:2] + 1j*data[1:taille:2];
      else:
        M = data;
      
      fileobj.close()

      V1 = np.reshape(M[0:nby*nbz], (nby, nbz));
      V2 = np.reshape(M[nby*nbz:nbx*nbz+nby*nbz], (nbx, nbz));
      V3 = np.reshape(M[nbz*(nbx+nby):nbz*(nbx+nby)+nbx*nby], (nbx, nby));
      if (raff != 1):
        V1 = Refine2D(ymin, ymax, zmin, zmax, V1, raff)
        V2 = Refine2D(xmin, xmax, zmin, zmax, V2, raff)
        V3 = Refine2D(xmin, xmax, ymin, ymax, V3, raff)

      X = np.linspace(xmin, xmax, (nbx-1)*raff+1);
      Y = np.linspace(ymin, ymax, (nby-1)*raff+1);
      Z = np.linspace(zmin, zmax, (nbz-1)*raff+1);
      XI = X
      YI = Y
      ZI = Z
      V = [np.transpose(V1),np.transpose(V2),np.transpose(V3)];
    elif (type_grid == 2):
      # SismoPlane
      offset = 0;
      xmin_ = np.zeros(nb_grids); xmax_ = np.zeros(nb_grids);
      ymin_ = np.zeros(nb_grids); ymax_ = np.zeros(nb_grids);
      zmin_ = np.zeros(nb_grids); zmax_ = np.zeros(nb_grids);
      centerx_ = np.zeros(nb_grids); centery_ = np.zeros(nb_grids);
      centerz_ = np.zeros(nb_grids);
      nbx_ = np.zeros(nb_grids, dtype='int'); nby_ = np.zeros(nb_grids, dtype='int');
      for i in range(nb_grids):
        xmin_[i] = np.fromfile(fileobj, prec, 1)[0]
        xmax_[i] = np.fromfile(fileobj, prec, 1)[0]
        ymin_[i] = np.fromfile(fileobj, prec, 1)[0]
        ymax_[i] = np.fromfile(fileobj, prec, 1)[0]
        zmin_[i] = np.fromfile(fileobj, prec, 1)[0]
        zmax_[i] = np.fromfile(fileobj, prec, 1)[0]
        centerx_[i] = np.fromfile(fileobj, prec, 1)[0]
        centery_[i] = np.fromfile(fileobj, prec, 1)[0]
        centerz_[i] = np.fromfile(fileobj, prec, 1)[0]
        nbx_[i] = np.fromfile(fileobj, 'i', 1)[0]
        nby_[i] = np.fromfile(fileobj, 'i', 1)[0]
        if (i < num_inst-1):
          offset = offset + nbx_[i]*nby_[i];
      
      origin = np.zeros(3); extA = np.zeros(3); extB = np.zeros(3)
      origin[0] = xmin_[num_inst-1]; extA[0] = xmax_[num_inst-1];
      origin[1] = ymin_[num_inst-1]; extA[1] = ymax_[num_inst-1];
      origin[2] = zmin_[num_inst-1]; extA[2] = zmax_[num_inst-1];
      extB[0] = centerx_[num_inst-1];
      extB[1] = centery_[num_inst-1]; extB[2] = centerz_[num_inst-1];
      nbx = nbx_[num_inst-1]; nby = nby_[num_inst-1];
      n = np.fromfile(fileobj, 'i', 1)[0]
      taille = cplx*nbx*nby;
      # we read values related to the required snapshot
      fileobj.seek(offset*prec_size, 1);
      data = np.fromfile(fileobj, prec, taille)
      if (cplx == 2):
        M = data[0:taille:2] + 1j*data[1:taille:2];
      else:
        M = data;
      
      fileobj.close()

      # we form the transpose of the original field (for display purposes)
      Tx = np.linspace(0, 1, (nbx-1)*raff+1);
      Ty = np.linspace(0, 1, (nby-1)*raff+1);
      [Tx, Ty] = np.meshgrid(Tx, Ty);
      U = np.reshape(M, (nby, nbx));
      if (raff != 1):
        U = Refine2D(0.0, 1.0, 0.0, 1.0, U, raff)
      
      V = U
      
      XI = (extA[0]-origin[0])*Tx + (extB[0]-origin[0])*Ty + origin[0];
      YI = (extA[1]-origin[1])*Tx + (extB[1]-origin[1])*Ty + origin[1];
      ZI = (extA[2]-origin[2])*Tx + (extB[2]-origin[2])*Ty + origin[2];
      
    elif (type_grid == 3):
      # SismoLine
      offset = 0;
      xmin_ = np.zeros(nb_grids); xmax_ = np.zeros(nb_grids);
      ymin_ = np.zeros(nb_grids); ymax_ = np.zeros(nb_grids);
      zmin_ = np.zeros(nb_grids); zmax_ = np.zeros(nb_grids);
      nbx_ = range(nb_grids)
      for i in range(nb_grids):
        xmin_[i] = np.fromfile(fileobj, prec, 1)[0]
        xmax_[i] = np.fromfile(fileobj, prec, 1)[0]
        ymin_[i] = np.fromfile(fileobj, prec, 1)[0]
        ymax_[i] = np.fromfile(fileobj, prec, 1)[0]
        zmin_[i] = np.fromfile(fileobj, prec, 1)[0]
        zmax_[i] = np.fromfile(fileobj, prec, 1)[0]
        nbx_[i] = np.fromfile(fileobj, 'i', 1)[0]
        if (i < num_inst-1):
          offset = offset + nbx_[i];
      
      xmin = xmin_[num_inst-1]; xmax = xmax_[num_inst-1];
      ymin = ymin_[num_inst-1]; ymax = ymax_[num_inst-1];
      zmin = zmin_[num_inst-1]; zmax = zmax_[num_inst-1];
      nbx = nbx_[num_inst-1];
      n = np.fromfile(fileobj, 'i', 1)[0]
      taille = cplx*nbx;
      # we read values related to the required snapshot
      fileobj.seek(offset*prec_size, 1);
      data = np.fromfile(fileobj, prec, taille)
      if (cplx == 2):
        M = data[0:taille:2] + 1j*data[1:taille:2];
      else:
        M = data;
      
      fileobj.close()
      
      T = np.linspace(0, 1, (nbx-1)*raff + 1); 
      V = M
      if (raff != 1):
        V = Refine1D(0.0, 1.0, M, raff)

      XI = np.linspace(xmin, xmax, nbx);
      YI = np.linspace(ymin, ymax, nbx);
      ZI = np.linspace(zmin, zmax, nbx);
      
    elif (type_grid == 4):
      # SismoPoint
      offset = 0;
      xmin_ = np.zeros(nb_grids);
      ymin_ = np.zeros(nb_grids);
      zmin_ = np.zeros(nb_grids);
      for i in range(nb_grids):
        xmin_[i] = np.fromfile(fileobj, prec, 1)[0]
        ymin_[i] = np.fromfile(fileobj, prec, 1)[0]
        zmin_[i] = np.fromfile(fileobj, prec, 1)[0]
        if (i < num_inst-1):
          offset = offset + 1;
      
      xmin = xmin_[num_inst-1];
      ymin = ymin_[num_inst-1];
      zmin = zmin_[num_inst-1];
      n = np.fromfile(fileobj, 'i', 1)[0]
      taille = cplx;
      # we read values related to the required snapshot
      fileobj.seek(offset*prec_size, 1);
      data = np.fromfile(fileobj, prec, taille)
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
      nbx_ = range(nb_grids)
      for i in range(nb_grids):
        nbx_[i] = np.fromfile(fileobj, 'i', 1)[0]
        coord = np.fromfile(fileobj, prec, 3*nbx_[i])
        if (i == num_inst-1):
          pts = coord
        if (i < num_inst-1):
          offset = offset + nbx_[i];

      nbx = int(nbx_[num_inst-1]);
      n = np.fromfile(fileobj, 'i', 1)[0]
      taille = cplx*nbx;
      # we read values related to the required snapshot
      fileobj.seek(offset*prec_size, 1);
      data = np.fromfile(fileobj, prec, taille)
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
      centerx_ = np.zeros(nb_grids)
      centery_ = np.zeros(nb_grids)
      centerz_ = np.zeros(nb_grids)
      normalex_ = np.zeros(nb_grids)
      normaley_ = np.zeros(nb_grids)
      normalez_ = np.zeros(nb_grids)
      radiusA_ = np.zeros(nb_grids)
      radiusB_ = np.zeros(nb_grids)
      nbx_ = range(nb_grids)
      for i in range(nb_grids):
        centerx_[i] = np.fromfile(fileobj, prec, 1)[0]
        centery_[i] = np.fromfile(fileobj, prec, 1)[0]
        centerz_[i] = np.fromfile(fileobj, prec, 1)[0]
        normalex_[i] = np.fromfile(fileobj, prec, 1)[0]
        normaley_[i] = np.fromfile(fileobj, prec, 1)[0]
        normalez_[i] = np.fromfile(fileobj, prec, 1)[0]
        radiusA_[i] = np.fromfile(fileobj, prec, 1)[0]
        radiusB_[i] = np.fromfile(fileobj, prec, 1)[0]
        nbx_[i] = np.fromfile(fileobj, 'i', 1)[0]
        if (i < num_inst-1):
          offset = offset + nbx_[i];

      centerx = centerx_[num_inst-1]; radiusA = radiusA_[num_inst-1];
      centery = centery_[num_inst-1]; radiusB = radiusB_[num_inst-1];
      centerz = centerz_[num_inst-1];
      normale = np.array([normalex_[num_inst-1], normaley_[num_inst-1], normalez_[num_inst-1]])
      nbx = nbx_[num_inst-1];
      n = np.fromfile(fileobj, 'i', 1)[0]
      taille = cplx*nbx;
      # we read values related to the required snapshot
      fileobj.seek(offset*prec_size, 1);
      data = np.fromfile(fileobj, prec, taille)
      if (cplx == 2):
        M = data[0:taille:2] + 1j*data[1:taille:2];
      else:
        M = data;
      
      fileobj.close()

      # finding the two vectors of plane, knowing the normale
      perm = np.array([0, 1, 2])
      if (abs(normale[0]) > max(normale[1],normale[2])):
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
      
      vec_u = np.zeros(3); vec_u[perm[2]] = -normale[perm[1]];
      vec_u[perm[1]] = normale[perm[2]]
      vec_v = np.cross(normale, vec_u)
      
      V = M
      if (raff != 1):
        V = Refine1D(0.0, 1.0, M, raff)

      teta = np.linspace(0, 2*pi, (nbx-1)*raff+2); teta = teta[0:(nbx-1)*raff+1];
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
    t = np.linspace(t[0], (taille-1)*dt + t[0], taille)
  
  E = np.zeros([taille, fin-debut]) + 1j*np.zeros([taille, fin-debut])    
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
  num = np.where(prod < 0)[0] + 1
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
    Vchap = np.fft.fft(V)

    # we keep only the signal that will contain the envelope
    if ((N0 == 0) and (N1 == 0)):
      num = GetLocalMaxima(Vchap, threshold*abs(Vchap).max())
      deb = (num[0] + num[1])//2
      Vchap[deb:] = 0
    else:
      Vchap[0:N0] = 0
      Vchap[N1:] = 0
    
    # returning back to real space
    Venv = np.fft.ifft(Vchap)
    return Venv
  elif (len(V.shape) == 2):
    N = V.shape[0]
    Venv = np.zeros(V.shape) + 1j*np.zeros(V.shape)
    for i in range(V.shape[1]):
      # Fourier transform is computed
      Vchap = np.fft.fft(V[:,i])

      # we keep only the signal that will contain the envelope
      if ((N0 == 0) and (N1 == 0)):
        num = GetLocalMaxima(Vchap, threshold*abs(Vchap).max())
        deb = (num[0] + num[1])//2
        Vchap[deb:] = 0
      else:
        Vchap[0:N0] = 0
        Vchap[N1:] = 0

      # returning back to real space
      Venv[:,i] = np.fft.ifft(Vchap)
    
    return Venv
  else:
    print("case not handled")
    
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
    x = np.reshape(args[2*i], np.size(args[2*i]))
    y = np.reshape(args[2*i+1], np.size(args[2*i+1]))
    droptol = threshold*max(abs(x))
    x = x*(abs(x) > droptol)*(abs(y) > droptol);
    y = y*(abs(x) > droptol)*(abs(y) > droptol);
    err += np.linalg.norm(x-y)**2
    norm_sol += np.linalg.norm(y)**2
    
  return np.sqrt(err/norm_sol)

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
  x = np.reshape(V,(np.size(V)))
  y = np.reshape(W,(np.size(V)))
  droptol = threshold*max(abs(x))
  x = x*(abs(x) > droptol)*(abs(y) > droptol);
  y = y*(abs(x) > droptol)*(abs(y) > droptol);
  diff = np.sort(abs(x-y))
  val = np.sort(abs(y))
  offset = min(np.where(val > 0)[0])
  N = np.size(x) - offset
  num = int(offset + ratio*N)
  err = diff[num]/val[num]
  return err

def erreurMode(V, mode):
  n = np.where(abs(mode)== abs(mode).max())[0]
  n1 = n[0]
  n2 = n[1]
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
  index = np.where(Delta<=0)[0]
  if(len(index) == 0):
    b = range(n)
    return b

  b = range(index[0]+1)
  for i in range(len(index)):
    iend = np.where(sig[index[i]+1:] > sig[index[i]-1])[0]
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
  V = np.loadtxt(nom_fichier)
  b = GetUniqueIndex(V[:,0])
  U = V[b,:]
  return U

def loadReverse(nom_fichier):
  V = np.loadtxt(nom_fichier)
  i = np.where(V[:,0] == V[:,0].min())[0]
  V[0:i+1, :] = np.flipud(V[0:i+1, :])
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
  taille = np.fromfile(fileobj, 'i', 1)[0]
  V = np.fromfile(fileobj, 'd', taille)
  fileobj.close()
  return V

def load1D_int(nom_fichier):
  """ loads a vector of integers (as produced by using Write function of Seldon vectors)
  Usage : V = load1D_int(nom_fichier)"""
  fileobj = open(nom_fichier, mode='rb')
  taille = np.fromfile(fileobj, 'i', 1)[0]
  V = np.fromfile(fileobj, 'i', taille)
  fileobj.close()
  return V

def load_full(nom_fichier, double_prec = True):
  """ loads a complex matrix (as produced by using Write function of Seldon dense matrices)
  Usage : A = load_full(nom_fichier)"""
  fileobj = open(nom_fichier, mode='rb')
  m = np.fromfile(fileobj, 'i', 1)[0]
  n = np.fromfile(fileobj, 'i', 1)[0]
  taille = m*n
  if (double_prec):
    Vz = np.fromfile(fileobj, 'd', 2*taille)
  else:
    Vz = np.fromfile(fileobj, 'f', 2*taille).astype('float64')
  
  V = Vz[0::2] + Vz[1::2]*1j;
  fileobj.close()
  A = np.reshape(V,(m, n))
  return A

def load_fullReal(nom_fichier, double_prec = True):
  """ loads a real matrix (as produced by using Write function of Seldon dense matrices)
  Usage : A = load_full(nom_fichier)"""
  fileobj = open(nom_fichier, mode='rb')
  m = np.fromfile(fileobj, 'i', 1)[0]
  n = np.fromfile(fileobj, 'i', 1)[0]
  taille = m*n
  if (double_prec):
    V = np.fromfile(fileobj, 'd', taille)
  else:
    V = np.fromfile(fileobj, 'f', taille).astype('float64')
  fileobj.close()
  A = np.reshape(V,(m, n))
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
  n = np.fromfile(fileobj, 'i', 1)[0]
  n = np.fromfile(fileobj, 'i', 1)[0]
  taille = (n+1)*n
  Vz = np.fromfile(fileobj, 'd', taille)
  V = Vz[0:taille:2] + Vz[1:taille:2]*1j;
  fileobj.close()
  A = np.zeros([n, n]) + 0j*np.zeros([n, n])
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
  A = np.loadtxt(nom_fichier)
  B = scipy.sparse.coo_matrix((A[:,2],(A[:,0]-1,A[:,1]-1)))
  return scipy.sparse.csr_matrix(B)
  
def loadComplexMat(nom_fichier):
  """ Reads a complex sparse matrix in coordinate format 
     Usage : A = loadComplexMat(nom_fichier) """
  fid = open(nom_fichier,'r')
  V = fid.readlines()
  fid.close()
  n = len(V)
  val = np.zeros(n) + 1j*np.zeros(n)
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
  m = np.fromfile(fileobj, 'i', 1)[0]
  n = np.fromfile(fileobj, 'i', 1)[0]
  k = np.fromfile(fileobj, 'i', 1)[0]
  taille = m*n*k
  if (double_prec):
    V = np.fromfile(fileobj, 'd', taille)
  else:
    V = np.fromfile(fileobj, 'f', taille).astype('float64')
  
  fileobj.close()
  A = np.reshape(V,(m, n, k))
  return A

def loadComplexVec(nom_fichier):
  """ Reads a complex vector in ascii format
     Usage : A = loadComplexVec(nom_fichier) """
  fid = open(nom_fichier,'r')
  V = fid.readlines()
  fid.close()
  n = len(V)
  val = np.zeros(n) + 1j*np.zeros(n)
  for i in range(n):
    ligne = V[i].strip();
    pos_comma = ligne.find(',')
    real_part = float(ligne[1:pos_comma])
    imag_part = float(ligne[pos_comma+1:len(ligne)-1])
    val[i] = real_part + 1j*imag_part
    
  return val

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
    v = np.array(v)/32767.0
    V = np.zeros([nframes, nchannel])
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
  type_data = np.fromfile(fid, 'i', 1)[0]
  prec = 'd'
  if (type_data == 0):
    prec = 'f';
    
  compt = 0
  while (compt < MAXI):
    try:
      n = np.fromfile(fid, 'i', 1)[0]
      Vloc = np.fromfile(fid, prec, n)
      if (compt == 0):
        L = np.zeros([MAXI, len(Vloc)])
      
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
  type_data = np.fromfile(fid, 'i', 1)[0]

  prec = 'd'; prec_size = 8;
  if ((type_data == 0) or (type_data == 2)):
    prec = 'f'
    prec_size = 4;
    
  n = np.fromfile(fid, 'i', 1)[0]
  Points = np.fromfile(fid, prec, n)
  Points = np.reshape(Points, [3,n/3])

  L = np.zeros([MAXI, n/3+1]);
  
  compt=0;
  while(compt<MAXI):
    try:
      n = np.fromfile(fid, 'i', 1)[0]
      Vloc = np.fromfile(fid, prec, n)
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

def ReadMeshData(fichier):
  fid = open(fichier, "r");
  type_data = np.fromfile(fid, 'i', 1)[0]
  print("type_data = ", type_data)
  nb_elt = np.fromfile(fid, 'i', 1)[0]
  offset = np.fromfile(fid, 'i', nb_elt+1)
  nu = [0]*nb_elt
  if (type_data == 0):
    for i in range(nb_elt):
      nu[i] = np.fromfile(fid, 'd', offset[i+1]-offset[i])
  else:
    for i in range(nb_elt):
      taille = 2*(offset[i+1]-offset[i])
      data = np.fromfile(fid, 'd', taille)
      nu[i] = data[0:taille:2] + 1j*data[1:taille:2]
  
  return nu

def WriteMeshData(nu, fichier):
  fid = open(fichier, "w")
  nb_elt = len(nu)
  if (nu[0].dtype == 'float64'):
    type_data = np.array([0, nb_elt])
  elif (nu[0].dtype == 'complex128'):
    type_data = np.array([1, nb_elt])
  
  type_data.tofile(fid)
  offset = np.array([0]*(nb_elt+1)).astype('int32')
  for i in range(nb_elt):
    offset[i+1] = offset[i] + len(nu[i])
    
  offset.tofile(fid)
  for i in range(nb_elt):
    nu[i].tofile(fid)
  
  fid.close()

def GetFourier(V, dt):
  """ Computes Fourier transform of V
  Usage : omega, Vhat = GetFourier(V, dt) 
  
          dt : time step between each value of V
          omega : pulsations for which Fourier transform is computed
          Vhat : Fourier transform for each point in omega  """
  
  coef_fft = dt / np.sqrt(2.0*pi)
  Vchap = coef_fft*np.fft.fft(V)
  omega = np.linspace(-pi/dt, pi/dt, len(V)+1)
  omega = omega[0:len(V)]
  return omega, np.fft.fftshift(Vchap)

def trace_fft(signal, F):
  """ Displays Fourier transform of a signal in decibels
   Usage : trace_fft(y, f)
   
   y : data to analyse
   f : sample frequency of the data  """
  
  fft_sig = np.fft.fft(signal);
  max_fft = max(abs(fft_sig));
  nfreq = len(signal)
  freq = np.linspace(-0.5*F, 0.5*F, len(signal)+1)
  np.plot(freq[0:nfreq], np.fft.fftshift(20.0*np.log10(abs(fft_sig)/max_fft)));
  np.xlim(0,0.4*F);
  np.xlabel('Frequency');
  np.ylabel('Fourier Transform : log(abs)')

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
    x = np.zeros(nb_vert)
    y = np.zeros(nb_vert)
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
    t = np.array(s[1:])
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
    m = max(m, np.size(A, 0));
    m = max(m, np.size(A, 1));
    
  
  A = np.zeros([m, m]) + 1j*np.zeros([m, m])
  for i in range(n):
    B = loadComplexMat(fic + '_P' + str(i) + '.dat').todense();
    if (assemble):
      A[0:np.size(B,0), 0:np.size(B,1)] += B;
    else:
      A[0:np.size(B,0), 0:np.size(B,1)] = B;

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

