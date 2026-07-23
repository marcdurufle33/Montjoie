#include "Quadrature/MontjoieQuadrature.hxx"

using namespace Montjoie;

template<class T>
void GetLobattoShapeFunctions(int r, Vector<UnivariatePolynomial<T> >& Phi,
			      Vector<UnivariatePolynomial<T> >& dPhi)
{
  T one(1), coef;
  
  Vector<T> xi, omega;
  ComputeGaussLobatto(xi, omega, r);
  
  // Lagrange polynomials with Gauss-Lobatto points
  Phi.Reallocate(r+1);
  UnivariatePolynomial<T> monome;
  monome.SetOrder(1);
  for (int i = 0; i <= r; i++)
    {
      Phi(i).SetOrder(0); Phi(i)(0) = one;
      for (int j = 0; j <= r; j++)
	if (i != j)
	  {
	    coef = one/(xi(i) - xi(j));
	    monome(0) = -xi(j)*coef;
	    monome(1) = coef;
	    Phi(i) = Phi(i)*monome;
	  }
    }

  //DISP(Phi);
  
  //for (int i = 0; i <= r; i++)
  //for (int j = 0; j <= r; j++)
  //DISP(Phi(i).Evaluate(xi(j)));  

  // derivatives of these functions
  dPhi.Reallocate(r+1);
  for (int i = 0; i <= r; i++)
    DerivatePolynomial(Phi(i), dPhi(i));
  
  //DISP(dPhi);
}

template<class T>
T ComputeSameElementMatrix(Vector<UnivariatePolynomial<T> >& Phi,
			   Vector<UnivariatePolynomial<T> >& dPhi,
			   const T& beta, int type, int n_ext, int n_int,
			   const UnivariatePolynomial<T>& jacob)
{
  T one(1), xi_i;
  
  // Gauss-Legendre rules
  int r = Phi.GetM()-1;
  Globatto<T> gauss_int, gauss_ext, gauss;
  gauss_int.ConstructQuadrature(n_int, gauss_int.QUADRATURE_GAUSS_SQUARED);
  //gauss_int.ConstructQuadrature(n_int, gauss_int.QUADRATURE_GAUSS);
  gauss_ext.ConstructQuadrature(n_ext);
  
  gauss.ConstructQuadrature(n_ext, gauss.QUADRATURE_GAUSS_SQUARED);
  VectReal_wp Points(2*n_ext+2), Weights(2*n_ext+2);
  for (int i = 0; i <= n_ext; i++)
    {
      xi_i = gauss.Points(i);
      Points(i) = xi_i/2;
      Weights(i) = gauss.Weights(i)/2;
      
      Points(i+n_ext+1) = (one-xi_i/2);
      Weights(i+n_ext+1) = gauss.Weights(i)/2;      
    }
  gauss_ext.AffectPoints(Points); gauss_ext.AffectWeights(Weights);
  n_ext = 2*n_ext+1;
  
  // valeur des fonctions de base sur les points de quadrature
  Matrix<T> GradPhiExt(r+1, n_ext+1), ValPhiExt(r+1, n_ext+1);
  for (int i = 0; i <= r; i++)
    for (int j = 0; j <= n_ext; j++)
      {
	ValPhiExt(i, j) = Phi(i).Evaluate(gauss_ext.Points(j));
	GradPhiExt(i, j) = dPhi(i).Evaluate(gauss_ext.Points(j));
      }

  Vector<Matrix<T> > GradPhiInt(n_ext+1), ValPhiInt(n_ext+1);
  Vector<Vector<T> > PointsInt(n_ext+1);
  Vector<Vector<T> > WeightsInt(n_ext+1);
  
  for (int i = 0; i <= n_ext; i++)
    {
      xi_i = gauss_ext.Points(i);
      if ((n_ext%2 == 0) && (i == n_ext/2))
	{
	  // cas d'un point au milieu, on ne met que deux sous-intervalles
	  PointsInt(i).Reallocate(2*(n_int+1));
	  WeightsInt(i).Reallocate(2*(n_int+1));
	  for (int j = 0; j <= n_int; j++)
	    {
	      PointsInt(i)(j) = (one-gauss_int.Points(j))/2;
	      PointsInt(i)(j+n_int+1) = (one+gauss_int.Points(j))/2;
	      WeightsInt(i)(j) = gauss_int.Weights(j)/2;
	      WeightsInt(i)(j+n_int+1) = gauss_int.Weights(j)/2;
	    }
	}
      else
	{
	  // cas avec juste deux intervalles [0, xi_i] et [xi_i, 1]
	  PointsInt(i).Reallocate(2*(n_int+1));
	  WeightsInt(i).Reallocate(2*(n_int+1));
	  for (int j = 0; j <= n_int; j++)
	    {
	      PointsInt(i)(j) = (one-gauss_int.Points(j))*xi_i;
	      PointsInt(i)(j+n_int+1) = gauss_int.Points(j)*(one-xi_i) + xi_i;
	      WeightsInt(i)(j) = gauss_int.Weights(j)*xi_i;
	      WeightsInt(i)(j+n_int+1) = gauss_int.Weights(j)*(one-xi_i);
	    }

	  // on met deux intervalles [0, xi_i] et [xi_i, 2xi_i] symetriques autour de xi_i
	  // et un intervalle [2 xi_i, 1]
	  /* PointsInt(i).Reallocate(3*(n_int+1));
	  WeightsInt(i).Reallocate(3*(n_int+1));
	  for (int j = 0; j <= n_int; j++)
	    {
	      PointsInt(i)(j) = (one-gauss_int.Points(j))*xi_i;
	      PointsInt(i)(j+n_int+1) = gauss_int.Points(j)*xi_i + xi_i;
	      PointsInt(i)(j+2*n_int+2) = gauss_int.Points(j)*(one-2*xi_i) + 2*xi_i;
	      WeightsInt(i)(j) = gauss_int.Weights(j)*xi_i;
	      WeightsInt(i)(j+n_int+1) = gauss_int.Weights(j)*xi_i;
	      WeightsInt(i)(j+2*n_int+2) = gauss_int.Weights(j)*(one-2*xi_i);
	      } */
	}
      
      // on evalue les fonctions sur les points de quadrature interieurs
      GradPhiInt(i).Reallocate(r+1, PointsInt(i).GetM());
      ValPhiInt(i).Reallocate(r+1, PointsInt(i).GetM());
      for (int k = 0; k <= r; k++)
	for (int j = 0; j < PointsInt(i).GetM(); j++)
	  {
	    GradPhiInt(i)(k, j) = dPhi(k).Evaluate(PointsInt(i)(j));
	    ValPhiInt(i)(k, j) = Phi(k).Evaluate(PointsInt(i)(j));
	  }
    }
  

  Matrix<T> mat_elem(r+1, r+1);

  // calcul de l'integrale singuliere
  // \int_0^1 \int_0^1 |s-t|^\beta f(s, t) ds dt

  // methode purement analytique
  // \int_0^1 (1-t)^\beta t^j dt = 1/ ( (beta+1+j) (beta+j) ... (beta+1) )
  Vector<T> IntElem(2*r+1+2*jacob.GetM()), InvBetaPlusN(2*r+3+2*jacob.GetM());
  for (int i = 0; i < InvBetaPlusN.GetM(); i++)
    InvBetaPlusN(i) = one/(beta+i+1);
  
  if (type == 0)
    {
      IntElem(0) = one/(beta+one);
      for (int i = 0; i < IntElem.GetM()-1; i++)
	IntElem(i+1) = IntElem(i)*(i+1)*InvBetaPlusN(i+1);
    }
  else
    {
      IntElem(0) = one/(beta+3);
      for (int i = 0; i < IntElem.GetM()-1; i++)
	IntElem(i+1) = IntElem(i)*(i+1)*InvBetaPlusN(i+3);    
    }
  
  Matrix<T> mat_ref(r+1, r+1);
  MultivariatePolynomial<T> PolI, PolJ, Q, PolJacob;
  UnivariatePolynomial<T> Pi, Pj;
  PolJacob.SetOrder(2, 2*jacob.GetOrder());
  for (int i = 0; i <= jacob.GetOrder(); i++)
    for (int j = 0; j <= jacob.GetOrder(); j++)
      PolJacob(i, j) = jacob(i)*jacob(j);
  
  for (int i = 0; i <= r; i++)
    for (int j = 0; j <= r; j++)
      {
	T val(0);
	if (type == 0)
	  {
	    if (jacob.GetOrder() == 0)
	      for (int k = 0; k < r; k++)
		for (int n = 0; n < r; n++)
		  val += dPhi(i)(k)*dPhi(j)(n)*InvBetaPlusN(k+n+1)*(IntElem(k)+IntElem(n));
	    else
	      {
		Pi = dPhi(i)*jacob; Pj = dPhi(j)*jacob;
		for (int k = 0; k <= Pi.GetOrder(); k++)
		  for (int n = 0; n <= Pj.GetOrder(); n++)
		    val += Pi(k)*Pj(n)*InvBetaPlusN(k+n+1)*(IntElem(k)+IntElem(n));
	      }
	  }
	else
	  {
	    PolI.SetOrder(2, r-1);
	    PolJ.SetOrder(2, r-1);
	    for (int n = 1; n <= r; n++)
	      for (int k = 0; k < n; k++)
		{
		  PolI(n-1-k, k) = Phi(i)(n);
		  PolJ(n-1-k, k) = Phi(j)(n);
		}
	    
	    if (jacob.GetOrder() == 0)
	      Q = PolI*PolJ;
	    else
	      Q = PolI*PolJ*PolJacob;
	    
	    for (int k = 0; k <= Q.GetOrder(); k++)
	      for (int n = 0; n <= Q.GetOrder()-k; n++)
		if (Q(k, n) != T(0))
		  val += Q(k, n)*InvBetaPlusN(k+n+3)*(IntElem(k)+IntElem(n));
	  }
	
	mat_ref(i, j) = val;
      }

  mat_elem.Fill(0);
  
  // methode purement numerique
  for (int k = 0; k <= n_ext; k++)
    {
      xi_i = gauss_ext.Points(k);
      T val_int(0);
      // avec Gauss-Squared
      for (int n = 0; n < PointsInt(k).GetM(); n++)
	{
	  val_int = gauss_ext.Weights(k)*WeightsInt(k)(n)*pow(abs(PointsInt(k)(n)-xi_i), beta);
	  if (jacob.GetOrder() > 0)
	    val_int *= jacob.Evaluate(xi_i)*jacob.Evaluate(PointsInt(k)(n));
	  
	  if (type == 0)
	    {
	      for (int i = 0; i <= r; i++)
		for (int j = 0; j <= r; j++)
		  mat_elem(i, j) += val_int*GradPhiInt(k)(i, n)*GradPhiExt(j, k);
	    }
	  else
	    {
	      for (int i = 0; i <= r; i++)
		for (int j = 0; j <= r; j++)
		  mat_elem(i, j) += val_int*(ValPhiInt(k)(i, n) - ValPhiExt(i, k))
		    *(ValPhiInt(k)(j, n) - ValPhiExt(j, k));
	    }
	}
    }
  
  //DISP(mat_ref);
  //DISP(mat_elem);
  
  Add(-one, mat_ref, mat_elem);
  cout << "Maximal error : " << MaxAbs(mat_elem) << endl;
  return MaxAbs(mat_elem);
}

template<class T>
T ComputeJointElementMatrix(Vector<UnivariatePolynomial<T> >& Phi,
			    Vector<UnivariatePolynomial<T> >& dPhi,
			    const T& beta, int type, int n_ext, int n_int,
			    const TinyVector<T, 2>& u, const TinyVector<T, 2>& v)
{
    // Gauss-Legendre rules
  int r = Phi.GetM()-1;
  int n_ref = 2*max(n_ext, n_int);
  Globatto<T> gauss_int, gauss_ext, gauss_ref, gauss_int_ext, gauss;
  gauss_int.ConstructQuadrature(n_int, gauss_int.QUADRATURE_GAUSS_SQUARED);
  //gauss_int.ConstructQuadrature(n_int, gauss_int.QUADRATURE_GAUSS);
  gauss_ext.ConstructQuadrature(n_ext);
  gauss_ref.ConstructQuadrature(n_ref);


  gauss.ConstructQuadrature(n_ext, gauss.QUADRATURE_GAUSS_SQUARED);
  VectReal_wp Points(2*n_ext+2), Weights(2*n_ext+2);
  for (int i = 0; i <= n_ext; i++)
    {
      T xi_i = gauss.Points(i), one(1);
      Points(i) = xi_i/2;
      Weights(i) = gauss.Weights(i)/2;
      
      Points(i+n_ext+1) = (one-xi_i/2);
      Weights(i+n_ext+1) = gauss.Weights(i)/2;      
    } 
  
  gauss_int_ext.AffectPoints(Points); gauss_int_ext.AffectWeights(Weights);
  
  //gauss_int_ext.ConstructQuadrature(n_ext);

  // valeur des fonctions de base sur les points de quadrature
  Matrix<T> GradPhiExt(r+1, n_ext+1), ValPhiExt(r+1, n_ext+1);
  for (int i = 0; i <= r; i++)
    for (int j = 0; j <= n_ext; j++)
      {
	ValPhiExt(i, j) = Phi(i).Evaluate(gauss_ext.Points(j));
	GradPhiExt(i, j) = dPhi(i).Evaluate(gauss_ext.Points(j));
      }

  int n_int_ext = gauss_int_ext.GetOrder();
  Matrix<T> GradPhiIntExt(r+1, n_int_ext+1), ValPhiIntExt(r+1, n_int_ext+1);
  for (int i = 0; i <= r; i++)
    for (int j = 0; j <= n_int_ext; j++)
      {
	ValPhiIntExt(i, j) = Phi(i).Evaluate(gauss_int_ext.Points(j));
	GradPhiIntExt(i, j) = dPhi(i).Evaluate(gauss_int_ext.Points(j));
      }

  Matrix<T> GradPhiInt(r+1, n_int+1), ValPhiInt(r+1, n_int+1);
  for (int i = 0; i <= r; i++)
    for (int j = 0; j <= n_int; j++)
      {
	ValPhiInt(i, j) = Phi(i).Evaluate(gauss_int.Points(j));
	GradPhiInt(i, j) = dPhi(i).Evaluate(gauss_int.Points(j));
      }

  T one(1), s, t;


  Matrix<T> mat_elem(r+1, r+1);

  // calcul de l'integrale singuliere
  // \int_0^1 \int_0^1 | (S + s u) - (S + t v) |^\beta f(s, t) ds dt
  Vector<T> IntElemS, IntElemT, InvBetaPlusN(2*r+3);
  for (int i = 0; i < InvBetaPlusN.GetM(); i++)
    InvBetaPlusN(i) = one/(beta+i+1);

  // methode semi-analytique avec deux fois plus de points
  // pour avoir une valeur de reference
  if (type == 0)
    {
      T val_int(0), Si;
      IntElemS.Reallocate(r); IntElemT.Reallocate(r);
      IntElemS.Fill(0); IntElemT.Fill(0);
      for (int k = 0; k <= n_ref; k++)
	{
	  s = gauss_ref.Points(k);
	  val_int = pow(sqrt(square(s*u(0)-v(0)) + square(s*u(1)-v(1))), beta)*gauss_ref.Weights(k);
	  Si = one;
	  for (int i = 0; i <= r-1; i++)
	    {
	      IntElemS(i) += val_int*Si;
	      Si *= s;
	    }
	  
	  t = s;
	  val_int = pow(sqrt(square(u(0)-t*v(0)) + square(u(1)-t*v(1))), beta)*gauss_ref.Weights(k);
	  Si = one;
	  for (int i = 0; i <= r-1; i++)
	    {
	      IntElemT(i) += val_int*Si;
	      Si *= t;
	    }
	}
    }
  else
    {
      IntElemS.Reallocate(2*r-1); IntElemT.Reallocate(2*r-1);
      T val_int(0), Si;
      IntElemS.Fill(0); IntElemT.Fill(0);
      for (int k = 0; k <= n_ref; k++)
	{
	  s = gauss_ref.Points(k);
	  val_int = pow(sqrt(square(s*u(0)-v(0)) + square(s*u(1)-v(1))), beta)*square(s-one)*gauss_ref.Weights(k);
	  Si = one;
	  for (int i = 0; i <= 2*r-2; i++)
	    {
	      IntElemS(i) += val_int*Si;
	      Si *= s;
	    }
	  
	  t = s;
	  val_int = pow(sqrt(square(u(0)-t*v(0)) + square(u(1)-t*v(1))), beta)*square(t-one)*gauss_ref.Weights(k);
	  Si = one;
	  for (int i = 0; i <= 2*r-2; i++)
	    {
	      IntElemT(i) += val_int*Si;
	      Si *= t;
	    }
	}
    }
  
  Matrix<T> mat_ref(r+1, r+1);
  MultivariatePolynomial<T> PolI, PolJ, Q;
  for (int i = 0; i <= r; i++)
    for (int j = 0; j <= r; j++)
      {
	T val(0);
	if (type == 0)
	  {
	    for (int k = 0; k < r; k++)
	      for (int n = 0; n < r; n++)
		val += dPhi(i)(k)*dPhi(j)(n)*InvBetaPlusN(k+n+1)*(IntElemS(k)+IntElemT(n));
	  }
	else
	  {
	    PolI.SetOrder(2, r-1);
	    PolJ.SetOrder(2, r-1);
	    for (int n = 1; n <= r; n++)
	      for (int k = 0; k < n; k++)
		{
		  PolI(n-1-k, k) = Phi(i)(n);
		  PolJ(n-1-k, k) = Phi(j)(n);
		}
	    
	    Q = PolI*PolJ;
	    for (int k = 0; k <= Q.GetOrder(); k++)
	      for (int n = 0; n <= Q.GetOrder()-k; n++)
		val += Q(k, n)*InvBetaPlusN(k+n+3)*(IntElemS(k)+IntElemT(n));
	  }
	
	mat_ref(i, j) = val;
      }

  // methode semi-analytique avec n_ext+1 points
  IntElemS.Fill(0); IntElemT.Fill(0);
  if (type == 0)
    {
      T val_int(0), Si;
      for (int k = 0; k <= n_ext; k++)
	{
	  s = gauss_ext.Points(k);
	  val_int = pow(sqrt(square(s*u(0)-v(0)) + square(s*u(1)-v(1))), beta)*gauss_ext.Weights(k);
	  Si = one;
	  for (int i = 0; i <= r-1; i++)
	    {
	      IntElemS(i) += val_int*Si;
	      Si *= s;
	    }
	  
	  t = s;
	  val_int = pow(sqrt(square(u(0)-t*v(0)) + square(u(1)-t*v(1))), beta)*gauss_ext.Weights(k);
	  Si = one;
	  for (int i = 0; i <= r-1; i++)
	    {
	      IntElemT(i) += val_int*Si;
	      Si *= t;
	    }
	}
    }
  else
    {
      T val_int(0), Si;
      for (int k = 0; k <= n_ext; k++)
	{
	  s = gauss_ext.Points(k);
	  val_int = pow(sqrt(square(s*u(0)-v(0)) + square(s*u(1)-v(1))), beta)*square(s-one)*gauss_ext.Weights(k);
	  Si = one;
	  for (int i = 0; i <= 2*r-2; i++)
	    {
	      IntElemS(i) += val_int*Si;
	      Si *= s;
	    }
	  
	  t = s;
	  val_int = pow(sqrt(square(u(0)-t*v(0)) + square(u(1)-t*v(1))), beta)*square(t-one)*gauss_ext.Weights(k);
	  Si = one;
	  for (int i = 0; i <= 2*r-2; i++)
	    {
	      IntElemT(i) += val_int*Si;
	      Si *= t;
	    }
	}
    }
  
  // methode semi-analytique
  Matrix<T> mat_semi(r+1, r+1);
  mat_semi.Fill(0);
  for (int i = 0; i <= r; i++)
    for (int j = 0; j <= r; j++)
      {
	T val(0);
	if (type == 0)
	  {
	    for (int k = 0; k < r; k++)
	      for (int n = 0; n < r; n++)
		val += dPhi(i)(k)*dPhi(j)(n)*InvBetaPlusN(k+n+1)*(IntElemS(k)+IntElemT(n));
	  }
	else
	  {
	    PolI.SetOrder(2, r-1);
	    PolJ.SetOrder(2, r-1);
	    for (int n = 1; n <= r; n++)
	      for (int k = 0; k < n; k++)
		{
		  PolI(n-1-k, k) = Phi(i)(n);
		  PolJ(n-1-k, k) = Phi(j)(n);
		}
	    
	    Q = PolI*PolJ;
	    for (int k = 0; k <= Q.GetOrder(); k++)
	      for (int n = 0; n <= Q.GetOrder()-k; n++)
		val += Q(k, n)*InvBetaPlusN(k+n+3)*(IntElemS(k)+IntElemT(n));
	  }
	
	mat_semi(i, j) = val;
      }

  // methode purement numerique
  mat_elem.Fill(0);
  for (int k = 0; k <= n_int_ext; k++)
    { 
      t = gauss_int_ext.Points(k);
      T val_int(0);
      for (int n = 0; n <= n_int; n++)
	{
	  s = gauss_int.Points(n);
	  val_int = gauss_int_ext.Weights(k)*gauss_int.Weights(n)*pow(sqrt(square(s*u(0)-t*v(0)) + square(s*u(1)-t*v(1))), beta);
	  
	  DISP(k); DISP(n); DISP(s); DISP(t); DISP(val_int);
	  if (type == 0)
	    {
	      for (int i = 0; i <= r; i++)
		for (int j = 0; j <= r; j++)
		  mat_elem(i, j) += val_int*GradPhiInt(i, n)*GradPhiIntExt(j, k);
	    }
	  else
	    {
	      for (int i = 0; i <= r; i++)
		for (int j = 0; j <= r; j++)
		  mat_elem(i, j) += val_int*(ValPhiInt(i, n) - ValPhiIntExt(i, k))*(ValPhiInt(j, n) - ValPhiIntExt(j, k));
	    }
	}
    }
  
  DISP(mat_ref);
  DISP(mat_semi);
  DISP(mat_elem);
  
  Add(-one, mat_ref, mat_semi);
  cout << "Maximal error ref : " << MaxAbs(mat_semi) << endl;

  Add(-one, mat_ref, mat_elem);
  cout << "Maximal error : " << MaxAbs(mat_elem) << endl;

  return MaxAbs(mat_elem);
}

Real_wp fac(int i)
{
  if (i == 0)
    return Real_wp(1);
  
  return i*fac(i-1);
}

int CheckTriangleQuadrature(int pmax, const VectR2& Points2D, const VectReal_wp& Weights2D)
{
  // testing integral of x^i y^j
  for (int p = 0; p <= pmax; p++)
    for (int i = 0; i <= p; i++)
      {
	int j = p-i;
        double val_exact = tgamma(i+1)*tgamma(j+1)/tgamma(i+j+3);
        double val = 0;
        for (int m = 0; m < Points2D.GetM(); m++)
          val += Weights2D(m)*pow(Points2D(m)(0), i)*pow(Points2D(m)(1), j);
        
        if (abs(val-val_exact)/val > 1e-12)
          {
	    return p-1;
            //cout << " Order " << i+j << " not exactly integrated for order p = " << p <<endl;
            //DISP(i); DISP(j); DISP(val); DISP(val_exact);
          }
      }  
  
  return pmax;
}

void GenerateMollerFormulaTriangle(int s)
{
  int d = 2*s+1, n = 2;
  VectR2 points; VectReal_wp poids;
  for (int i = 0; i <= s; i++)
    {
      Real_wp w = pow(2.0, -2.0*s)*pow(d+n-2*i, d) / (fac(i)*fac(d+n-i));
      if (i%2 == 1)
        w = -w;
      
      for (int beta1 = 0; beta1 <= s-i; beta1++)
        for (int beta2 = 0; beta2 <= s-i-beta1; beta2++)
          {
            Real_wp x = Real_wp(2*beta1+1) / (d + n - 2*i);
            Real_wp y = Real_wp(2*beta2+1) / (d + n - 2*i);
            points.PushBack(R2(x, y));
            poids.PushBack(w);
          }
    }

  DISP(points); DISP(poids);
  cout << "Number of points = " << points.GetM() << endl;
  int p = CheckTriangleQuadrature(d+1, points, poids);
  cout << "Degree = " << p << endl;

  TriangleQuadrature::ConstructQuadrature(d, points, poids, TriangleQuadrature::QUADRATURE_GAUSS);
  cout << "Number of gaussian points = " << points.GetM() << endl;
}

int main(int argc, char** argv)
{
  if (argc < 2)
    {
      cout << "Number of arguments too small" << endl;
      abort();
    }
  
  InitMontjoie(argc, argv);
  
  int r = atoi(argv[1]);
  GenerateMollerFormulaTriangle(r);

  return 0;
  
  //typedef long double T;
  typedef Real_wp T;
  DISP(GetPrecisionMachine(T(1)));
  
  Vector<UnivariatePolynomial<T> > Phi, dPhi;
  GetLobattoShapeFunctions(r, Phi, dPhi);
  
  UnivariatePolynomial<T> jacob(3);
  //jacob.FillRand(); Mlt(Real_wp(1e-9), jacob); DISP(jacob);
  jacob.SetOrder(0); jacob(0) = 1.0;
  Vector<int> array_n;
  //array_n.PushBack(1); array_n.PushBack(2); array_n.PushBack(3); array_n.PushBack(4);
  //array_n.PushBack(6); array_n.PushBack(8); array_n.PushBack(12); array_n.PushBack(16);
  //array_n.PushBack(20); array_n.PushBack(22); array_n.PushBack(24);
  array_n.PushBack(26);
  //array_n.PushBack(30); array_n.PushBack(40); array_n.PushBack(50); array_n.PushBack(70);
  //array_n.PushBack(100);
  //array_n.PushBack(200);   array_n.PushBack(300);   array_n.PushBack(500);   //array_n.PushBack(1000);
  T err(1), err_prev; int n_int, n_ext;
  //TinyVector<T, 2> u (1.0, -0.2), v(-1.0, 0.19);
  TinyVector<T, 2> u(1.46446609406726, -3.53553390593274), v(-3.53553390593274, 1.46446609406726);
  for (int i = 0; i < array_n.GetM(); i++)
  //for (int i = 0; i < 3*r+100; i+=2)
    {
      n_int = r+11; n_ext = r+array_n(i);
      //n_int = array_n(i);
      //n_ext = array_n(i);
      //n_ext = array_n(i);
      //n_int *= r; n_ext *= r;
      DISP(n_int); DISP(n_ext);
      
      err_prev = err;
      //err = ComputeSameElementMatrix(Phi, dPhi, T(-0.5), 0, n_ext, n_int, jacob);
      err = ComputeJointElementMatrix(Phi, dPhi, T(0.5), 0, n_ext, n_int, u, v);
      if (i > 0)
      {
      cout << "Order = " << (log(err) - log(err_prev))/(log(T(n_int)) - log(T(r*array_n(i-1)))) << endl;
      }						      
    }


  return FinalizeMontjoie();
}
