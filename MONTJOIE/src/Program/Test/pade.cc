#include "Montjoie.hxx"

using namespace Montjoie;

Real_wp fac(int m)
{
  return tgamma(Real_wp(m+1));
}

// computes roots of a polynomial P, real parts of roots are contained in R
// imaginary roots are contained in Rimag
Real_wp GetRoots(const UnivariatePolynomial<Real_wp>& P, VectReal_wp& R, VectReal_wp& Rimag)
{
  int first_coef = P.GetOrder();
  for (int i = 0; i < P.GetM(); i++)
    if (abs(P(i)) > 1e3*epsilon_machine)
      {
	first_coef = i;
	break;
      }
  
  int last_coef = 0;
  Real_wp one(1), zero(0), coef;
  for (int i = P.GetM()-1; i >= 0; i--)
    if (P(i) != zero)
      {
	last_coef = i;
	break;
      }
  
  // companion matrix
  int N = last_coef - first_coef;
  Matrix<Real_wp> A, V;
  A.Reallocate(N, N);
  A.Fill(zero);
  for (int i = 0; i < N-1; i++)
    A(i+1, i) = one;
  
  coef = one/P(last_coef);
  for (int i = first_coef; i < last_coef; i++)
    A(i-first_coef, N-1) = -coef*P(i);
  
  DISP(A);
  // roots are equal to eigenvalues of the companion matrix
  GetEigenvaluesEigenvectors(A, R, Rimag, V);
  
  return P(first_coef);
}

int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);

  if (argc != 2)
    {
      cout << "Entrez m" << endl;
      abort();
    }

  int m = to_num<int>(argv[1]);
  
  // calcul des racines de P
  UnivariatePolynomial<Real_wp> P;
  P.SetOrder(m);
  for (int i = 0; i <= m; i ++)
    {
      Real_wp coef = fac(m)*fac(2*m-i) / (fac(2*m) * fac(i) * fac(m-i));
      P(i) = coef;
    }

  DISP(P);
  VectReal_wp R, Rimag;
  GetRoots(P, R, Rimag);

  VectReal_wp CoefC(2*m+1);
  for (int k = 0; k <= 2*m; k++)
    CoefC(k) = Real_wp(1)/(pow(Real_wp(2), k-1) * fac(k));

  DISP(CoefC);
  // Calcul des ci : on prend les points de Gauss-Lobatto
  VectReal_wp ci, weights;
  if (m == 1)
    {
      ci.Reallocate(1);
      ci(0) = Real_wp(0.5);
    }
  else
    {
      //ComputeGaussLobatto(ci, weights, 2*m-2);
      ComputeGaussLegendre(ci, weights, m-1);
    }
  
  ci = ci - Real_wp(0.5);
  DISP(ci);
  
  // Matrice de Vandermonde pour les ci
  // les factorielles sont mises dans le second membre
  Matrix<Real_wp> VDM_all(2*m-1, ci.GetM());  
  for (int i = 0; i < 2*m-1; i++)
    for (int j = 0; j < ci.GetM(); j++)
      VDM_all(i, j) = pow(ci(j), i);

  Matrix<Real_wp> VDM(VDM_all);
  VDM.Resize(ci.GetM(), ci.GetM());
  DISP(VDM); DISP(VDM_all);
  
  // on la factorise
  Vector<int> pivot;
  GetLU(VDM, pivot);  

  // boucle sur les puissances de L
  for (int powL = 0; powL < m; powL++)
    {
      VectReal_wp vec_f(2*m);
      vec_f.Zero();
      for (int i = 0; i <= m; i += 2)
	{
	  Real_wp rho_IM = P(i);
	  for (int k = 1; k < 2*m-i; k += 2)
	    {
	      int j = k+i - powL;
	      if ((j >= 1) && (j <= k))
		{
		  Real_wp ck = CoefC(k);
		  vec_f(j-1) += rho_IM*ck;
		}
	    }
	}
      
      for (int i = 1; i <= m; i += 2)
	{
	  Real_wp rho_IM = P(i);
	  for (int k = 0; k < 2*m-i; k += 2)
	    {
	      int j = k+i - powL;
	      if ((j >= 1) && (j <= k))
		{		  
		  Real_wp ck = CoefC(k);
		  vec_f(j-1) -= rho_IM * ck;
		}
	    }
	}

      DISP(powL); DISP(vec_f);

      // coefficients avec forme alternative
      int r = powL+1;
      int jmin = 2;
      if (r%2 == 1)
	jmin = 1;

      vec_f.Zero();
      for (int j = jmin; j <= 2*m-r; j += 2)
	{
	  //DISP(j);
	  for (int i = 0; i <= min(m, r-1); i += 2)
	    {
	      //DISP(i); DISP(P(i)); DISP(r+j-i-1); DISP(CoefC(r+j-i-1));
	      vec_f(j-1) += P(i)*CoefC(r+j-i-1); //DISP(vec_f(j-1));
	    }

	  for (int i = 1; i <= min(m, r-1); i += 2)
	    {
	      //DISP(i); DISP(P(i)); DISP(r+j-i-1); DISP(CoefC(r+j-i-1));
	      vec_f(j-1) -= P(i)*CoefC(r+j-i-1);	  //DISP(vec_f(j-1));
	    }
	}

      DISP(vec_f);

      // forme alternative finale
      vec_f.Zero();
      for (int j = 0; j <= 2*m-r; j++)
	{
	  if (j%2 != r%2)
	    for (int i = 0; i <= min(m, r-1); i++)
	      vec_f(j) += pow(-1.0, i)*P(i)*CoefC(r+j-i);
	}

      DISP(vec_f);
      
      VectReal_wp rhs(ci.GetM());
      rhs.Zero();
      for (int i = 0; i < ci.GetM(); i++)
	rhs(i) = fac(i)*vec_f(i);
      
      DISP(rhs);
      SolveLU(VDM, pivot, rhs);
      
      cout << "Pour powL = " << powL << endl;
      cout << "wi = " << rhs << endl;
      cout << endl << endl;

      VectReal_wp vec_w(2*m-1); vec_w.Zero();
      Mlt(VDM_all, rhs, vec_w);

      for (int i = 0; i < 2*m-1; i++)
	vec_w(i) /= fac(i);
      
      DISP(vec_w);

    }

  PadeScheme_Iterator<Real_wp> pade;
  pade.SetOrder(2*m);

  DISP(pade.GetWeights());

  return 0;
}
