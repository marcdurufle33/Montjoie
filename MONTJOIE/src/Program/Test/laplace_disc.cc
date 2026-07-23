#include "Montjoie.hxx"

using namespace Montjoie;

VectComplex_wp ComputePhi(const VectReal_wp& eps, const VectReal_wp& sigma, const Real_wp& omega, const VectReal_wp& radius)
{
  int N = radius.GetM() - 1;
  VectComplex_wp rhs(2*N);
  Matrix<Complex_wp> mat(2*N, 2*N);
  mat.Zero(); rhs.Zero();

  mat(0, 0) = log(radius(0));
  mat(0, 1) = 1.0; rhs(0) = 1.0;
  for (int i = 0; i < N-1; i++)
    {
      mat(2*i+1, 2*i) = log(radius(i+1));
      mat(2*i+1, 2*i+1) = 1.0;
      mat(2*i+1, 2*i+2) = -log(radius(i+1));
      mat(2*i+1, 2*i+3) = -1.0;
      mat(2*i+2, 2*i) = eps(i) + Iwp*sigma(i)/omega;
      mat(2*i+2, 2*i+2) = -eps(i+1) - Iwp*sigma(i+1)/omega;
    }

  mat(2*N-1, 2*N-2) = log(radius(N));
  mat(2*N-1, 2*N-1) = 1.0;

  Vector<int> pivot;
  GetLU(mat, pivot);

  SolveLU(mat, pivot, rhs);
  return rhs;
}

int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);

  int N = 3;
  VectReal_wp eps(N), sigma(N), radius(N); Real_wp omega = 2e-7*2*pi_wp;
  eps(0) = 1.0; eps(1) = 1.0; eps(2) = 81.0;
  sigma(0) = 0.0; sigma(1) = 4.346e8; sigma(2) = 1883.65;
  radius(0) = 0.0175; radius(1) = 0.092; radius(2) = 0.1; radius(3) = 20.0;
  
  VectComplex_wp coef = ComputePhi(eps, sigma, omega, radius);

  DISP(coef);
  
  
  return FinalizeMontjoie();
}
