#include "Share/MontjoieCommon.hxx"

using namespace Montjoie;

int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);
  
  VarRandomGenerator var;
  int n = 245, nx, ny, nz;
  var.RoundToSquare(n, 1.5, nx, ny);
  DISP(nx); DISP(ny); DISP(n);

  n = 2456;
  var.RoundToSquare(n, 1.5, 1.8, nx, ny, nz);
  DISP(nx); DISP(ny); DISP(nz); DISP(n);
  
  DISP(var.GetRand());
  
  VectReal_wp x;
  var.GenerateRandomNumbers(20, x);
  DISP(x);
  
  IVect permut;
  var.GenerateRandomPermutation(20, permut);
  DISP(permut);
  var.ApplyRandomPermutation(x);
  DISP(x);
  
  var.SetRandomGenerator(var.GSL_MT19937);
  var.GenerateRandomNumbers(20, x);
  DISP(x);
  var.GenerateRandomPermutation(20, permut);
  DISP(permut);
  
  IVect prime_n;
  FindPrimeNumbers(prime_n, 20);
  DISP(prime_n);

  int p = 2*2*2*3*3*11*17*17;
  DISP(p);
  IVect prime_decomp;
  FindPrimeDecomposition(p, prime_decomp, prime_n);
  DISP(prime_decomp);

  FindPrimeNumbers(prime_n, 5);
  p = 37;
  FindPrimeFactorization(p, prime_decomp, prime_n);
  DISP(p); DISP(prime_decomp);
  
  n = 20;
  FindTwoFactors(n, nx, ny);
  DISP(nx); DISP(ny);
  
  int nb_points = FindClosestPow2357(975);
  DISP(nb_points);
  FindPrimeDecomposition(nb_points, prime_decomp, prime_n);
  DISP(prime_decomp);

  IVect decomp;
  DecomposeContinuedFraction(pi_wp, Real_wp(1e-6), decomp); 
  
  // then you can compute p and q, such that p/q is the fraction approximating pi
  int q;
  GetNumeratorDenominator(decomp, p, q);
  DISP(p); DISP(q);

  return FinalizeMontjoie();
}
