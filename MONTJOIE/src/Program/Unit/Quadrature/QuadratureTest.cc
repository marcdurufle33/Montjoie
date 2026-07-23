#define MONTJOIE_WITH_TWO_DIM
#define MONTJOIE_WITH_THREE_DIM

#include "Quadrature/MontjoieQuadrature.hxx"

using namespace Montjoie;

Real_wp threshold;

bool EqualQuadrature(const Globatto<Real_wp>& lob, const VectReal_wp& points,
		     const VectReal_wp& weights)
{
  if ((lob.GetOrder() != points.GetM()-1) || (lob.GetGeometryOrder() != points.GetM()-1))
    {
      cout << "GetOrder/GetGeometryOrder incorrect" << endl;
      abort();
    }
  
  if (lob.GetNbPointsQuad() != points.GetM())
    {
      cout << "GetNbPointsQuad incorrect" << endl;
      abort();
    }
  
  for (int i = 0; i <= lob.GetOrder(); i++)
    if (isnan(points(i)) || isnan(weights(i)) || isnan(lob.Points(i)) || isnan(lob.Weights(i)) ||
	(abs(points(i)-lob.Points(i)) > 0.1*threshold) || (abs(weights(i)-lob.Weights(i)) > 0.1*threshold) ||
	(abs(points(i)-lob.Points()(i)) > 0.1*threshold) || (abs(weights(i)-lob.Weights()(i)) > 0.1*threshold) )
      {
	cout << "Points/Weights incorrect" << endl;
	return false;
	//abort();
      }
  
  return true;
}

int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);
  
  threshold = 1e4*epsilon_machine;
  
  /* glob_chrono.Start(VirtualTimer::ALL);
  
  int r = atoi(argv[1]);
  Globatto<Real_wp> lob_t;
  lob_t.ConstructQuadrature(r);
  if (r < 10)
    {
      cout << "Points = " << lob_t.Points() << endl;
      cout << "Weights = " << lob_t.Weights() << endl;
    }

  glob_chrono.Stop(VirtualTimer::ALL);
  cout << "Temps ecoule = " << glob_chrono.GetSeconds(VirtualTimer::ALL) << endl;

  exit(0); */
  
  // evaluation of recurrence relation for orthogonal polynomials
  // for example : Legendre polynomials 
  Matrix<Real_wp> ab;
  GetJacobiPolynomial(ab, 5, Real_wp(0), Real_wp(0));
  cout<<"P0 = 1"<<endl;
  cout<<"P1 = x - "<<ab(0,0)<<endl;
  for (int k = 1; k < ab.GetM(); k++) 
    cout<<"P"<<k+1<<" = (x - "<<ab(k,0)<<") P"<<k<<" - "<<ab(k, 1)<<" P"<<k-1<<endl;
  
  // then you can evaluate those polynomials
  Vector<Real_wp> Pn, dPn;
  EvaluateJacobiPolynomial(ab, 5, Real_wp(0.23), Pn);
  cout<<"Legendre polynomials at x = 0.23 "<<endl<<Pn<<endl;
  
  // Legendre polynomials with direct evaluation
  Vector<UnivariatePolynomial<Real_wp> > Ln;
  ComputeLegendrePolynome(0, 7, Ln);
  for (int i = 0; i < 6; i++)
    {
      Real_wp Pn_ref = Ln(i).Evaluate(Real_wp(0.23)) / Ln(i)(Ln(i).GetOrder());
      if (isnan(Pn_ref) || isnan(Pn(i)) || (abs(Pn_ref - Pn(i)) > threshold))
	{
	  cout << " GetJacobiPolynomial/EvaluateJacobiPolynomial incorrect" << endl;
	  abort();
	}
    }

  EvaluateJacobiPolynomial(ab, 5, Real_wp(0.23), Pn, dPn);
  for (int i = 0; i < 6; i++)
    {
      UnivariatePolynomial<Real_wp> dLn;
      DerivatePolynomial(Ln(i), dLn);
      Real_wp dPn_ref = dLn.Evaluate(Real_wp(0.23)) / Ln(i)(Ln(i).GetOrder());
      if (isnan(dPn_ref) || isnan(dPn(i)) || (abs(dPn_ref - dPn(i)) > threshold))
	{
	  cout << "EvaluateJacobiPolynomial incorrect" << endl;
	  abort();
	}
    }
  
  // you can find roots of last polynomial (P6)
  Vector<Real_wp> roots, w;
  SolveGauss(roots, w, 5, ab);
  cout<<"Roots of Legendre Polynomials "<<roots<<endl;
  
  for (int i = 0; i < roots.GetM(); i++)
    if (abs(Ln(5).Evaluate(2.0*roots(i)-1.0)) > threshold)
      {
	cout << "SolveGauss incorrect" << endl;
	abort();
      }
  
  // you can construct Gauss points and weights directly :
  Vector<Real_wp> points, weights;
  for (int r = 1; r <= 10; r++)
    {
      ComputeGaussLegendre(points, weights, r);
      Vector<Real_wp> points_gauss(points), weights_gauss(weights);
      //cout<<"Gauss points for r = 3 : "<<endl<<points<<endl;
      //cout<<"Gauss weights for r = 3 : "<<endl<<weights<<endl;
      // checking the accuracy
      int order_int = 2*r+2;
      for (int p = 0; p <= 2*r+2; p++)
	{
	  // evaluating integral of x^p on [0,1]
	  Real_wp val = 0;
	  for (int i = 0; i < points.GetM(); i++)
	    val += weights(i)*pow(points(i), p);
	  
	  Real_wp val_exact = Real_wp(1)/(p+1);
	  if (abs(val - val_exact) > 0.1*val_exact*threshold)
	    {
	      order_int = p-1;
	      //cout<<"Order "<<p<<" is not exactly integrated "<<endl;
	      break;
	    }
	}
      
      if (order_int != (2*r+1))
	{
	  cout << "ComputeGaussLegendre incorrect" << endl;
	  abort();
	}
  
      // same stuff for Gauss-Lobatto (two extremities are quadrature points)
      ComputeGaussLobatto(points, weights, r);
      Vector<Real_wp> points_lobatto(points), weights_lobatto(weights);
      //cout<<"Gauss-Lobatto points for r = 3 : "<<endl<<points<<endl;
      //cout<<"Gauss-Lobatto weights for r = 3 : "<<endl<<weights<<endl;
      // checking the accuracy
      order_int = 2*r+2;
      for (int p = 0; p <= 2*r+2; p++)
	{
	  // evaluating integral of x^p on [0,1]
	  Real_wp val = 0;
	  for (int i = 0; i < points.GetM(); i++)
	    val += weights(i)*pow(points(i), p);
	  
	  Real_wp val_exact = Real_wp(1)/(p+1);
	  if (abs(val - val_exact) > 0.1*val_exact*threshold)
	    {
	      //cout<<"Order "<<p<<" is not exactly integrated "<<endl;
	      order_int = p-1;
	      break;
	    }
	}
      
      if (order_int != (2*r-1))
	{
	  cout << "ComputeGaussLobatto incorrect" << endl;
	  abort();
	}
      
      // Gauss-Blended formulas
      ComputeGaussBlendedFormulas(r+1, points, weights, 0.0);
      for (int i = 0; i < points.GetM(); i++)
	if ((abs(points(i) - points_gauss(i)) > 0.1*threshold)
	    || isnan(points(i)) || isnan(points_gauss(i))
	    || (abs(weights(i) - weights_gauss(i)) > 0.1*threshold))
	  {
	    cout << "ComputeGaussBlendedFormulas incorrect" << endl;
	    abort();
	  }

      ComputeGaussBlendedFormulas(r+1, points, weights, 1.0);
      for (int i = 0; i < points.GetM(); i++)
	if ((abs(points(i) - points_lobatto(i)) > 0.1*threshold)
	    || isnan(points(i)) || isnan(points_lobatto(i))
	    || (abs(weights(i) - weights_lobatto(i)) > 0.1*threshold))
	  {
	    cout << "ComputeGaussBlendedFormulas incorrect" << endl;
	    abort();
	  }

      ComputeGaussBlendedFormulas(r+1, points, weights, 7.0/8);
      order_int = 2*r+2;
      for (int p = 0; p <= 2*r+2; p++)
	{
	  // evaluating integral of x^p on [0,1]
	  Real_wp val = 0;
	  for (int i = 0; i < points.GetM(); i++)
	    val += weights(i)*pow(points(i), p);
	  
	  Real_wp val_exact = Real_wp(1)/(p+1);
	  if (abs(val - val_exact) > 0.1*val_exact*threshold)
	    {
	      //cout<<"Order "<<p<<" is not exactly integrated "<<endl;
	      order_int = p-1;
	      break;
	    }
	}
      
      if (order_int != (2*r-1))
	{
	  cout << "ComputeGaussBlendedFormulas incorrect" << endl;
	  abort();
	}
      
      // same stuff for Gauss-Jacobi
      ComputeGaussJacobi(points, weights, r, Real_wp(2), Real_wp(1));
      //cout<<"Gauss-Jacobi points for r = 3 (alpha = 2, beta = 1) : "<<endl<<points<<endl;
      //cout<<"Gauss-Jacobi weights for r = 3 (alpha = 2, beta = 1) : "<<endl<<weights<<endl;
      // checking the accuracy
      order_int = 2*r+2;
      for (int p = 0; p <= 2*r+2; p++)
	{
	  // evaluating \int_0^1 (1-x)^2 x x^p  = 2! p+1! / (p+4)! = 2 / ((p+2)(p+3)(p+4)) 
	  Real_wp val = 0;
	  for (int i = 0; i < points.GetM(); i++)
	    val += weights(i)*pow(points(i), p);
	  
	  Real_wp val_exact = Real_wp(2)/((p+2)*(p+3)*(p+4));
	  if (abs(val - val_exact) > 0.1*val_exact*threshold)
	    {
	      order_int = p-1;
	      //cout<<"Order "<<p<<" is not exactly integrated "<<endl;
	      break;
	    }
	}

      if (order_int != (2*r+1))
	{
	  cout << "ComputeGaussJacobi incorrect" << endl;
	  abort();
	}
      
      // same stuff for Gauss-Radau-Jacobi (an extremity is a quadrature point)
      // true -> 1.0 is the extremity where the quadrature point is placed
      ComputeGaussRadauJacobi(points, weights, r, Real_wp(2), Real_wp(1), true);
      //cout<<"Gauss-Radau-Jacobi points for r = 3 (alpha = 2, beta = 1) : "<<endl<<points<<endl;
      //cout<<"Gauss-Radau-Jacobi weights for r = 3 (alpha = 2, beta = 1) : "<<endl<<weights<<endl;
      // checking the accuracy
      order_int = 2*r+2;
      for (int p = 0; p <= 2*r+2; p++)
	{
	  // evaluating \int_0^1 (1-x)^2 x x^p  = 2! p+1! / (p+4)! = 2 / ((p+2)(p+3)(p+4)) 
	  Real_wp val = 0;
	  for (int i = 0; i < points.GetM(); i++)
	    val += weights(i)*pow(points(i), p);
	  
	  Real_wp val_exact = Real_wp(2)/((p+2)*(p+3)*(p+4));
	  if (abs(val - val_exact) > 0.1*val_exact*threshold)
	    {
	      order_int = p-1;
	      //cout<<"Order "<<p<<" is not exactly integrated "<<endl;
	      break;
	    }
	}

      if ((order_int != (2*r)) || (abs(points(r)-1.0) > 0.1*threshold))
	{
	  cout << "ComputeGaussRadauJacobi incorrect" << endl;
	  abort();
	}
      
      // same stuff for Gauss-Lobatto-Jacobi (two extremities are quadrature points)
      ComputeGaussLobattoJacobi(points, weights, r, Real_wp(2), Real_wp(1));
      //cout<<"Gauss-Lobatto-Jacobi points for r = 3 (alpha = 2, beta = 1) : "<<endl<<points<<endl;
      //cout<<"Gauss-Lobatto-Jacobi weights for r = 3 (alpha = 2, beta = 1) : "<<endl<<weights<<endl;
      // checking the accuracy
      order_int = 2*r+2;
      for (int p = 0; p <= 2*r+2; p++)
	{
	  // evaluating \int_0^1 (1-x)^2 x x^p  = 2! p+1! / (p+4)! = 2 / ((p+2)(p+3)(p+4)) 
	  Real_wp val = 0;
	  for (int i = 0; i < points.GetM(); i++)
	    val += weights(i)*pow(points(i), p);
	  
	  Real_wp val_exact = Real_wp(2)/((p+2)*(p+3)*(p+4));
	  if (abs(val - val_exact) > 0.1*val_exact*threshold)
	    {
	      order_int = p-1;
	      //cout<<"Order "<<p<<" is not exactly integrated "<<endl;
	      break;
	    }
	}
      
      if ((order_int != (2*r-1)) || (abs(points(r)-1.0) > 0.1*threshold) || (abs(points(0)) > 0.1*threshold))
	{
	  cout << "ComputeGaussLobattoJacobi incorrect" << endl;
	  abort();
	}

    }
 
  // testing Gauss-Chebyshev
  ComputeGaussChebyshev(points, weights, 6);
  Real_wp one(1);
  int order_int = 14;
  for (int p = 0; p <= 14; p++)
    {
      Real_wp val_ref = 0;
      switch(p)
	{
	case 0 : val_ref = pi_wp; break;
	case 2 : val_ref = pi_wp/2; break;
	case 4 : val_ref = 3*pi_wp/8; break;
	case 6 : val_ref = 5*pi_wp/16; break;
	case 8 : val_ref = 35*pi_wp/128; break;
	case 10 : val_ref = 63*pi_wp/256; break;
	case 12 : val_ref = 231*pi_wp/1024; break;
	case 14 : val_ref = 429*pi_wp/2048; break;
	}
      
      Real_wp val = 0;
      for (int i = 0; i < points.GetM(); i++)
	{
	  Real_wp x = 2*points(i)-one;
	  val += 2*pow(x, p)*weights(i);
	}
      
      if (abs(val - val_ref) > 0.1*threshold)
	{
	  order_int = p-1;
	  //cout << "Order " << p << " is not exactly integrated " << endl;
	  break;
	}
    }
  
  if (order_int != 13)
    {
      cout << "ComputeGaussChebyshev incorrect" <<endl;
      abort();
    }
  
  // testing Gauss-Logarithmic
  Real_wp a = -0.5;
  ComputeGaussLogarithmic(points, weights, 9, a);
  order_int = 20;
  for (int p = 0; p <= 20; p++)
    {
      Real_wp val_ref = one/square(p+1+a);
      
      Real_wp val = 0;
      for (int i = 0; i < points.GetM(); i++)
	val += pow(points(i), p)*weights(i);
      
      if (abs(val - val_ref) > 0.1*val_ref*threshold)
	{
	  //DISP(p); DISP(val-val_ref);
	  order_int = p-1;
	  //cout << "Order " << p << " is not exactly integrated " << endl;
	  break;
	}
    }

  if (order_int != 19)
    {
      DISP(order_int);
      cout << "ComputeGaussLogarithmic incorrect" <<endl;
      abort();
    }

  // testing Gauss-Laguerre
  a = Real_wp(-0.5);
  ComputeGaussLaguerre(points, weights, 4, a);
  order_int = 10;
  for (int p = 0; p <= 10; p++)
    {
      Real_wp val_ref = tgamma(p+a+1);
      
      Real_wp val = 0;
      for (int i = 0; i < points.GetM(); i++)
	val += pow(points(i), p)*weights(i);
      
      if (abs(val - val_ref)/val_ref > 0.1*val_ref*threshold)
	{
	  order_int = p-1;
	  //DISP(val-val_ref);
	  //cout << "Order " << p << " is not exactly integrated " << endl;
	  break;
	}
    }

  if (order_int != 9)
    {
      DISP(order_int);
      cout << "ComputeGaussLaguerre incorrect" <<endl;
      abort();
    }
  
  // testing Gauss-Hermite
  a = Real_wp(2);
  ComputeGaussHermite(points, weights, 4, a);
  order_int = 10;
  for (int p = 0; p <= 10; p++)
    {
      Real_wp val_ref(0);
      if ((p+toInteger(a))%2 == 0)
	{
	  val_ref = sqrt(pi_wp);      
	  for (int k = 2; k <= p+toInteger(a); k += 2)
	    val_ref *= Real_wp(k-1)/2;
	}
      
      Real_wp val = 0;
      for (int i = 0; i < points.GetM(); i++)
	val += pow(points(i), p)*weights(i);
      
      if (abs(val - val_ref) > 0.1*threshold)
	{
	  order_int = p-1;
	  //DISP(val-val_ref);
	  //cout << "Order " << p << " is not exactly integrated " << endl;
	  break;
	}
    }

  if (order_int != 9)
    {
      DISP(order_int);
      cout << "ComputeGaussHermite incorrect" <<endl;
      abort();
    }
  
  // testing Globatto in 1-D
  Globatto<Real_wp> lob;
  Globatto<Real_wp>::blending_default = 0.4;
  
  for (int r = 1; r <= 10; r++)
    {
      lob.ConstructQuadrature(r);
      ComputeGaussLegendre(points, weights, r);
      if (!EqualQuadrature(lob, points, weights))
	{
	  cout << "ConstructQuadrature(gauss) incorrect" << endl;
	  abort();
	}
  
      ComputeGaussLobatto(points, weights, r);
      lob.ConstructQuadrature(r, lob.QUADRATURE_LOBATTO);
      if (!EqualQuadrature(lob, points, weights))
	{
	  cout << "ConstructQuadrature(lobatto) incorrect" << endl;
	  abort();
	}

      Real_wp a = 0.5, b = 2.5;
      ComputeGaussRadauJacobi(points, weights, r, Real_wp(0), Real_wp(0), false);
      lob.ConstructQuadrature(r, lob.QUADRATURE_RADAU);
      if (!EqualQuadrature(lob, points, weights))
	{
	  cout << "ConstructQuadrature(radau) incorrect" << endl;
	  abort();
	}

      ComputeGaussJacobi(points, weights, r, a, b);
      lob.ConstructQuadrature(r, lob.QUADRATURE_JACOBI, a, b);
      if (!EqualQuadrature(lob, points, weights))
	{
	  cout << "ConstructQuadrature(jacobi) incorrect" << endl;
	  abort();
	}

      ComputeGaussRadauJacobi(points, weights, r, a, b, false);
      lob.ConstructQuadrature(r, lob.QUADRATURE_RADAU_JACOBI, a, b);
      if (!EqualQuadrature(lob, points, weights))
	{
	  cout << "ConstructQuadrature(radau_jacobi) incorrect" << endl;
	  abort();
	}

      ComputeGaussLobattoJacobi(points, weights, r, a, b);
      lob.ConstructQuadrature(r, lob.QUADRATURE_LOBATTO_JACOBI, a, b);
      if (!EqualQuadrature(lob, points, weights))
	{
	  cout << "ConstructQuadrature(lobatto_jacobi) incorrect" << endl;
	  abort();
	}
      
      lob.ConstructQuadrature(r, lob.QUADRATURE_GAUSS_SQUARED);
      order_int = 2*r+2;
      for (int p = 0; p <= 2*r+2; p++)
	{
	  Real_wp val = 0;
	  for (int i = 0; i < lob.GetNbPointsQuad(); i++)
	    val += lob.Weights(i)*pow(lob.Points(i), Real_wp(p-1)/2);
	  
	  Real_wp val_ref = Real_wp(2)/(p+1);
	  if (isnan(val) || isnan(val_ref) || abs(val-val_ref) > 0.1*val_ref*threshold)
	    {
	      order_int = p-1;
	      break;
	    }	  
	}
      
      if (order_int != 2*r+1)
	{
	  cout << "ConstructQuadrature(gauss_squared) incorrect" << endl;
	  abort();
	}

      ComputeGaussBlendedFormulas(r+1, points, weights, Globatto<Real_wp>::blending_default);
      lob.ConstructQuadrature(r, lob.QUADRATURE_GAUSS_BLENDED);
      if (!EqualQuadrature(lob, points, weights))
	{
	  cout << "ConstructQuadrature(gauss_blended) incorrect" << endl;
	  abort();
	}
    }      

  // you can select your own interpolation points
  points.Reallocate(7);
  points(0) = Real_wp(0.01); points(1) = 0.12; points(2) = 0.238; 
  points(3) = 0.289; points(4) = 0.421; points(5) = 0.718; points(6) = Real_wp(0.992);
  
  Sort(points);
  
  weights.Reallocate(points.GetM()); weights.FillRand(); Mlt(1e-9, weights);
  lob.AffectPoints(points); lob.AffectWeights(weights);
  for (int i = 0; i < points.GetM(); i++)
    if (isnan(points(i)) || isnan(lob.Points(i)) || (abs(points(i)-lob.Points(i)) > 0.1*threshold))
      {
	cout << "AffectPoints incorrect" << endl;
	abort();
      }

  for (int i = 0; i < weights.GetM(); i++)
    if (isnan(weights(i)) || isnan(lob.Weights(i)) || (abs(weights(i)-lob.Weights(i)) > 0.1*threshold))
      {
	cout << "AffectWeights incorrect" << endl;
	abort();
      }
  
  // checking basis functions
  Vector<UnivariatePolynomial<Real_wp> > phi(lob.GetOrder()+1), dphi(lob.GetOrder()+1);
  for (int i = 0; i <= lob.GetOrder(); i++)
    {
      phi(i).SetOrder(0);
      phi(i)(0) = Real_wp(1);
      
      UnivariatePolynomial<Real_wp> factor;
      factor.SetOrder(1);
      for (int j = 0; j <= lob.GetOrder(); j++)
	if (j != i)
	  {
	    factor(1) = 1.0/(lob.Points(i)-lob.Points(j));
	    factor(0) = -factor(1)*lob.Points(j);
	    phi(i) *= factor;
	  }
      
      DerivatePolynomial(phi(i), dphi(i));
    }
  
  Linspace(0.0, 1.0, 100, points);
  for (int i = 0; i < points.GetM(); i++)
    {
      lob.ComputeValuesPhiRef(points(i), Pn);
      for (int j = 0; j <= lob.GetOrder(); j++)
	{
	  Real_wp val = lob.EvaluatePhi(j, points(i));
	  Real_wp dval = lob.EvaluatePhiGrad(j, points(i));
	  Real_wp val_ref = phi(j).Evaluate(points(i));
	  Real_wp dval_ref = dphi(j).Evaluate(points(i));
	  if (isnan(val) || isnan(val_ref) || abs(val-val_ref) > 100*threshold)
	    {
	      DISP(val); DISP(val_ref); DISP(i); DISP(j); DISP(points(i));
	      cout << "EvaluatePhi incorrect" << endl;
	      abort();
	    }
	  
	  val = Pn(j);
	  if (isnan(val) || isnan(val_ref) || abs(val-val_ref) > 100*threshold)
	    {
	      DISP(val); DISP(val_ref); DISP(i); DISP(j); DISP(points(i));
	      cout << "ComputeValuesPhiRef incorrect" << endl;
	      abort();
	    }

	  if (isnan(dval) || isnan(dval_ref) || abs(dval-dval_ref) > 100*threshold)
	    {
	      DISP(dval); DISP(dval_ref); DISP(i); DISP(j); DISP(points(i));
	      cout << "EvaluatePhiGrad incorrect" << endl;
	      abort();
	    }
	}

    }
  
  // testing ComputeGradPhi
  lob.ComputeGradPhi();
  
  for (int i = 0; i <= lob.GetOrder(); i++)
    for (int j = 0; j <= lob.GetOrder(); j++)
      {
	Real_wp val_ref = dphi(i).Evaluate(lob.Points(j));
	if (isnan(lob.GradPhi(i, j)) || isnan(val_ref) 
	    || isnan(lob.GradPhi()(i, j)) || (abs(lob.GradPhi(i, j) - val_ref) > 100*threshold) ||
	    (abs(lob.GradPhi()(i, j) - val_ref) > 100*threshold) )
	  {
	    DISP(i); DISP(j); DISP(val_ref); DISP(lob.GradPhi(i, j)); DISP(lob.GradPhi()(i, j));
	    cout << "ComputeGradPhi incorrect" << endl;
	    abort();
	  }
      }
  
  // testing SubdivGlobatto
  int nb_subdiv = 5, order = 4;
  lob.ConstructQuadrature(order, lob.QUADRATURE_LOBATTO);
  Vector<Globatto<Real_wp> > vec_phi(nb_subdiv);
  Vector<Real_wp> all_points(1);
  all_points(0) = Real_wp(0);
  for (int i = 0; i < nb_subdiv; i++)
    {
      Real_wp x0 = Real_wp(i)/nb_subdiv;
      Real_wp xf = Real_wp(i+1)/nb_subdiv;
      w.Reallocate(order+1);
      for (int j = 0; j <= order; j++)
	{
	  w(j) = x0 + (xf-x0)*lob.Points(j);
	  if (j > 0)
	    all_points.PushBack(w(j));
	}
      
      vec_phi(i).AffectPoints(w);
    }
  
  SubdivGlobatto lob_s;
  lob_s.Init(false, nb_subdiv, order);
  if (lob_s.GetOrder() != nb_subdiv*order)
    {
      cout << "GetOrder incorrect" << endl;
      abort();
    }

  if (lob_s.GetGeometryOrder() != nb_subdiv*order)
    {
      cout << "GetGeometryOrder incorrect" << endl;
      abort();
    }
  
  for (int i = 0; i <= lob_s.GetOrder(); i++)
    if (isnan(all_points(i)) || isnan(lob_s.Points(i)) || (abs(all_points(i) - lob_s.Points(i)) > threshold))
      {
	cout << "Points incorrect" << endl;
	abort();
      }
  
  for (int j = 0; j < points.GetM(); j++)
    {
      int n = 0;
      for (int k = 0; k < nb_subdiv; k++)
	if ((points(j) >= (all_points(k*order)-threshold))
	    && (points(j) <= (all_points((k+1)*order)+threshold)))
	  {
	    n = k;
	    break;
	  }
      
      for (int i = 0; i <= lob_s.GetOrder(); i++)
	{
	  Real_wp val = lob_s.EvaluatePhi(i, points(j));
	  Real_wp val_ref = 0;
	  int iloc = -1;
	  if ((i >= n*order) && (i <= (n+1)*order))
	    iloc = i-n*order;
	  
	  if (iloc >= 0)
	    val_ref = vec_phi(n).EvaluatePhi(iloc, points(j));
	  
	  if (isnan(val) || isnan(val_ref) || (abs(val-val_ref) > threshold))
	    {
	      cout << "EvaluatePhi incorrect" << endl;
	      abort();
	    }
	}
    }
  
  // 2-D quadrature points
  /*VectR2 Points2D; VectReal_wp Weights2D;
  TriangleQuadrature::ConstructQuadrature(3, Points2D, Weights2D,
					  TriangleQuadrature::QUADRATURE_GAUSS);
  
  cout<<"2-D points " << endl << Points2D << endl;
  cout<<"2-D weights " << endl << Weights2D << endl;

  VectR3 Points3D; VectReal_wp Weights3D;
  TetrahedronQuadrature::ConstructQuadrature(3, Points3D, Weights3D,
					     TetrahedronQuadrature::QUADRATURE_GAUSS);
  
  cout << "3-D points " << endl << Points3D << endl;
  cout << "3-D weights " << endl << Weights3D << endl; */
  
  cout << "All tests passed successfully" << endl;
  
  return FinalizeMontjoie();
}
