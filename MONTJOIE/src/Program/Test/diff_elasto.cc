#include "Montjoie.hxx"

using namespace Montjoie;

int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);

  int Nx = 100, Ny = 100;
  Matrix<Real_wp, General, ArrayRowSparse> A;
  VectReal_wp b;
  
  A.Reallocate(2*Nx*Ny, 2*Nx*Ny);
  b.Reallocate(2*Nx*Ny);
  for (int i = 0; i < Nx; i++)
    for (int j = 0; j < Ny; j++)
      {
      }

  

  return FinalizeMontjoie();
}
