#include "Solver/MontjoieSolver.hxx"

using namespace Montjoie;

class RungeKuttaPolynom
{
public :
  UnivariatePolynomial<Real_wp> ESin, ECos;
  
  void SetOrder(int r, int N = 0)
  {
    UnivariatePolynomial<Real_wp> monome;
    Real_wp one(1), coef;
    UnivariatePolynomial<Real_wp> x, pone;
    x.SetOrder(1); pone.SetOrder(0);
    pone(0) = one;
    x(1) = one;
    if (N == 0)
      {
	// pas de stabilisation
	ECos = pone; ESin = x;
	coef = -one/2;
	monome = x*x;
	for (int i = 2; i <= r; i += 2)
	  {
	    Add(coef, monome,  ECos);
	    coef *= -one/((i+1)*(i+2));
	    monome = monome*x*x;
	  }
	
	monome = x*x*x;
	coef = -one/6;
	for (int i = 3; i <= r; i += 2)
	  {
	    Add(coef, monome,  ESin);
	    coef *= -one/((i+1)*(i+2));
	    monome *= x*x;      
	  }
      }
    else
      {
	UnivariatePolynomial<Real_wp> Tn, Tnm1, Tnm2, Tnp1, X;
	Real_wp beta = sqrt(Real_wp((r-1)*(r-1)-1));
	X = x; X(1) /= beta;
	Tnm2 = pone;
	Tnm1 = X;
	Tn = 2.0*X*Tnm1 - Tnm2;
	for (int n = 3; n <= r; n++)
	  {
	    Tnp1 = 2.0*X*Tn - Tnm1;
	    Tnm2 = Tnm1;
	    Tnm1 = Tn;
	    Tn = Tnp1;
	  }
	
	Tn *= 0.5*(r-2)/sqrt(beta*beta+1);
	Tnm2 *= -0.5*r/sqrt(beta*beta+1);
	Tnm1 *= beta/sqrt(beta*beta+1);
	DISP(Tnm1); DISP(Tn); DISP(Tnm2);
	
	ESin.SetOrder(r-1);
	ECos.SetOrder(r);
	for (int k = 1; k <= r-1; k += 2)
	  ESin(k) = Tnm1(k);
	
	ECos(r) = Tn(r);
	for (int k = 0; k <= r-2; k += 2)
	  ECos(k) = Tn(k)+Tnm2(k);
	
	if (r%4 == 2)
	  ECos *= -one;
	else
	  ESin *= -one;
	
	DISP(ECos); DISP(ESin);
      }
  }

};

Real_wp GetRoots(const UnivariatePolynomial<Real_wp>& P, VectReal_wp& R, VectReal_wp& Rimag)
{
  int first_coef = P.GetOrder();
  for (int i = 0; i < P.GetM(); i++)
    if (abs(P(i)) > epsilon_machine)
      {
	first_coef = i;
	break;
      }
  
  //cout << "First non-null coefficient = " << P(first_coef)
  //<< " occuring for x^" << first_coef << endl;
  
  int last_coef = 0;
  Real_wp one(1), zero(0), coef;
  for (int i = P.GetM()-1; i >= 0; i--)
    if (P(i) != zero)
      {
	last_coef = i;
	break;
      }
  
  int N = last_coef - first_coef;
  if (N <= 0)
    return Real_wp(0);
  
  Matrix<Real_wp> A, V;
  A.Reallocate(N, N);
  A.Fill(zero);
  for (int i = 0; i < N-1; i++)
    A(i+1, i) = one;
  
  coef = one/P(last_coef);
  for (int i = first_coef; i < last_coef; i++)
    A(i-first_coef, N-1) = -coef*P(i);
  
  GetEigenvaluesEigenvectors(A, R, Rimag, V);
  
  return P(first_coef);
}


Real_wp GetCFL(const UnivariatePolynomial<Real_wp>& P)
{
  VectReal_wp roots, root_imag;
  Real_wp first_coef = GetRoots(P, roots, root_imag);
  if (first_coef > 0)
    return Real_wp(0);
  
  Real_wp cfl = 1e10; DISP(roots);
  for (int i = 0; i < roots.GetM(); i++)
    if ((root_imag(i) == Real_wp(0)) && (roots(i) < -epsilon_machine) && (abs(roots(i)) < cfl))
      cfl = abs(roots(i));
  
  return cfl;
}

int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);
 
  if (argc < 2)
    {
      cout << "Entrez l'ordre" << endl;
      return -1;
    }
  
  int K = atoi(argv[1]);
  
  ofstream file_out("cfl.dat");
  file_out.precision(15);
  // calcul de la cfl a tout ordre
  for (int r = 1; r <= K; r++)
    {
      //if ((r%4 ==1)||(r%4==2))
      if ((r%2 ==1))
	{
	  // unstable scheme
	}
      else
	{
	  RungeKuttaPolynom P;
	  P.SetOrder(r, 4);
	  UnivariatePolynomial<Real_wp> Module, pone, ValReal;
	  pone.SetOrder(0); pone(0) = Real_wp(1) + Real_wp(1e-15);
	  Module = P.ECos*P.ECos + P.ESin*P.ESin - pone;
	  ValReal = P.ECos + P.ESin - pone;
	  for (int i = 1; i < ValReal.GetM(); i++)
	    if (ValReal(i) < 0)
	      ValReal(i) = -ValReal(i);
	  
	  Real_wp cfl = GetCFL(Module);
	  Real_wp cfl_real = GetCFL(ValReal);
	  	  
	  cout << "Order = " << r << endl;	  
	  cout << "CFL = " << cfl << endl;
	  file_out << "r = " << r << " , CFL = " << cfl << "CFL on real axis = " << cfl_real << endl;
	  file_out << "ECos = poly1d([" << P.ECos(P.ECos.GetOrder());
	  for (int i = P.ECos.GetOrder()-1; i >= 0; i--)
	    file_out << ", " << P.ECos(i);
	  
	  file_out << "])" << endl;
	  file_out << "ESin = poly1d([" << P.ESin(P.ESin.GetOrder());
	  for (int i = P.ESin.GetOrder()-1; i >= 0; i--)
	    file_out << ", " << P.ESin(i);
	  
	  file_out << "])" << endl;
	  
	}
    }
    
  
  // calcul de la cfl avec un parametre
  Real_wp a = 0, b = 0.01, one(1);
  int N = 200;
  RungeKuttaPolynom P;
  P.SetOrder(K);
  UnivariatePolynomial<Real_wp> Module, pone, monome, x;
  x.SetOrder(1); x(1) = one;
  pone.SetOrder(0); pone(0) = Real_wp(1)+ Real_wp(1e-20);
  monome = Pow(x, K+1);
  for (int i = 0; i < N; i++)
    {
      Real_wp lambda = Real_wp(i)/(N-1);
      Real_wp coef = (one-lambda)*a + lambda*b;
      
      Module = square(P.ECos) + square(P.ESin+coef*monome) - pone;
      Real_wp cfl = GetCFL(Module);
      file_out << coef << " " << cfl << endl;
    }

  // calcul de la cfl avec deux parametres
  /* Real_wp a1 = 1e-10, b1 = 0.001, a2 = 1e-10, b2 = 0.0002,  one(1);
  int Nx = 301, Ny = 301;
  RungeKuttaPolynom P;
  P.SetOrder(K);
  UnivariatePolynomial<Real_wp> Module, pone, monome1, monome2, x;
  x.SetOrder(1); x(1) = one;
  pone.SetOrder(0); pone(0) = Real_wp(1)+ Real_wp(1e-20);
  monome1 = Pow(x, K+1);
  monome2 = Pow(x, K+2);
  for (int i = 0; i < Nx; i++)
    for (int j = 0; j < Ny; j++)
      {
	Real_wp lambda = Real_wp(i)/(Nx-1);
	Real_wp mu = Real_wp(j)/(Ny-1);
	Real_wp coef1 = (one-lambda)*a1 + lambda*b1;
	Real_wp coef2 = (one-mu)*a2 + mu*b2;
      
	Module = square(P.ECos-coef2*monome2) + square(P.ESin+coef1*monome1) - pone;
	Real_wp cfl = GetCFL(Module);
	file_out << coef1 << " " << coef2 << " " << cfl << endl;
      }
  */
  
  file_out.close();
  
  return FinalizeMontjoie();
}
