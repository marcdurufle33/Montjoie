#include "Share/MontjoieCommon.hxx"

using namespace Montjoie;

Real_wp threshold = 1e-12;

Real_wp fac(int l)
{
  Real_wp res(1);
  for (int i = 2; i <= l; i++)
    res *= i;

  return res;
}

template<class PolynomialTest>
void CheckPolynomial(PolynomialTest& Ptest)
{
  PolynomialTest T;
  
  // on teste le constructeur et destructeur
  if (T.GetOrder() != 0)
    {
      cout << "Constructeur par defaut incorrect" << endl;
      abort();
    }
  
  int n = 5;
  PolynomialTest P(n);
  
  if (P.GetOrder() != n)
    {
      cout << "Constructeur avec le degre du polynome incorrect" << endl;
      abort();
    }
  
  // on teste l'operateur d'acces ()
  P(0) = 2.0; P(1) = 0.5; P(2) = -2.3; P(3) = 2.7; P(4) = -0.6; P(5) = 1.8;
  
  if ( (P(0) != 2.0) || (P(1) != 0.5) || (P(2) != -2.3) || (P(3) != 2.7) || (P(4) != -0.6) || (P(5) != 1.8))
    {
      cout << "operateur d'acces incorrect" << endl;
      abort();
    }
  
  // on teste l'operateur <<
  cout << "P = " << P << endl;
  
  // on teste l'operateur =
  PolynomialTest Q = P;
  
  for (int i = 0; i <= n; i++)
    if (Q(i) != P(i))
      {
	cout << "operateur = incorrect" << endl;
	abort();
      }
  
  // on teste SetOrder
  PolynomialTest Ps;
  Ps.SetOrder(n);

  if (Ps.GetOrder() != n)
    {
      cout << "SetOrder incorrect" << endl;
      abort();
    }

  Ps(0) = 2.0; Ps(1) = 0.5; Ps(2) = -2.3; Ps(3) = 2.7; Ps(4) = -0.6; Ps(5) = 1.8;
  
  if ( (Ps(0) != 2.0) || (Ps(1) != 0.5) || (Ps(2) != -2.3) || (Ps(3) != 2.7)
       || (Ps(4) != -0.6) || (Ps(5) != 1.8))
    {
      cout << "SetOrder incorrect" << endl;
      abort();
    }
  
  // on teste ResizeOrder
  Ps.ResizeOrder(n+4);
  
  Ps(6) = 1.2; Ps(7) = -3.3; Ps(8) = 4.7; Ps(9) = -5.1;
  if ( (Ps(0) != 2.0) || (Ps(1) != 0.5) || (Ps(2) != -2.3) || (Ps(3) != 2.7)
       || (Ps(4) != -0.6) || (Ps(5) != 1.8) || (Ps(6) != 1.2) || (Ps(7) != -3.3) 
       || (Ps(8) != 4.7) || (Ps(9) != -5.1) )
    {
      cout << "ResizeOrder incorrect" << endl;
      abort();
    }

  // on teste la fonction Evaluate
  Seldon::Vector<double> x;
  Linspace(0.0, 1.0, 100, x);
  for (int i = 0; i < x.GetM(); i++)
    {
      double val = P.Evaluate(x(i));
      double val_ex = 0.0;
      for (int j = 0; j <= n; j++)
	val_ex += pow(x(i), j)*P(j);
      
      if (abs(val-val_ex) > threshold)
	{
	  cout << "Fonction Evaluate incorrecte" << endl;
	  abort();
	}
    }
  
  // operateur += pour des polynomes
  for (int i = 0; i <= n; i++)
    Q(i) = double(rand())/RAND_MAX;
  
  PolynomialTest Plow(n-2), Phigh(n+3);
  for (int i = 0; i <= Plow.GetOrder(); i++)
    Plow(i) = double(rand())/RAND_MAX;

  for (int i = 0; i <= Phigh.GetOrder(); i++)
    Phigh(i) = double(rand())/RAND_MAX;
  
  P = Q;
  P += Plow;
  
  for (int i = 0; i < x.GetM(); i++)
    {
      double val = P.Evaluate(x(i));
      double res = Q.Evaluate(x(i)) + Plow.Evaluate(x(i));
      if ((abs(val -res) > threshold) || isnan(val) || isnan(res))
	{
	  cout << "operateur += incorrect" << endl;
	  abort();
	}
    }

  P = Q;
  P += Phigh;
  
  for (int i = 0; i < x.GetM(); i++)
    {
      double val = P.Evaluate(x(i));
      double res = Q.Evaluate(x(i)) + Phigh.Evaluate(x(i));
      if ((abs(val -res) > threshold) || isnan(val) || isnan(res))
	{
	  cout << "operateur += incorrect" << endl;
	  abort();
	}
    }

  // operateur -= pour des polynomes
  P = Q;
  P -= Plow;
  
  for (int i = 0; i < x.GetM(); i++)
    {
      double val = P.Evaluate(x(i));
      double res = Q.Evaluate(x(i)) - Plow.Evaluate(x(i));
      if ((abs(val -res) > threshold) || isnan(val) || isnan(res))
	{
	  cout << "operateur -= incorrect" << endl;
	  abort();
	}
    }

  P = Q;
  P -= Phigh;
  
  for (int i = 0; i < x.GetM(); i++)
    {
      double val = P.Evaluate(x(i));
      double res = Q.Evaluate(x(i)) - Phigh.Evaluate(x(i));
      if ((abs(val -res) > threshold) || isnan(val) || isnan(res))
	{
	  cout << "operateur -= incorrect" << endl;
	  abort();
	}
    }
  
  // operateur *= pour des polynomes
  P = Q;
  P *= Plow;
  
  for (int i = 0; i < x.GetM(); i++)
    {
      double val = P.Evaluate(x(i));
      double res = Q.Evaluate(x(i)) * Plow.Evaluate(x(i));
      if ((abs(val -res) > threshold) || isnan(val) || isnan(res))
	{
	  cout << "operateur *= incorrect" << endl;
	  abort();
	}
    }

  P = Q;
  P *= Phigh;
  
  for (int i = 0; i < x.GetM(); i++)
    {
      double val = P.Evaluate(x(i));
      double res = Q.Evaluate(x(i)) * Phigh.Evaluate(x(i));
      if ((abs(val -res) > threshold) || isnan(val) || isnan(res))
	{
	  cout << "operateur *= incorrect" << endl;
	  abort();
	}
    }
  
  P = Q;
  
  // on teste l'operateur +=
  Q += 0.5;
  if (abs(Q(0)-P(0)-0.5) > threshold)
    {
      cout << "operateur += incorrect" << endl;
      abort();
    }
  
  for (int i = 1; i <= n; i++)
    if (Q(i) != P(i))
      {
	cout << "operateur += incorrect" << endl;
	abort();
      }
  
  // on teste l'operateur -=
  Q -= 0.2;
  if (abs(Q(0)-P(0)-0.3) > threshold)
    {
      cout << "operateur -= incorrect" << endl;
      abort();
    }
  
  for (int i = 1; i <= n; i++)
    if (Q(i) != P(i))
      {
	cout << "operateur -= incorrect" << endl;
	abort();
      }
  
  // on teste l'operateur *=
  Q(0) = P(0);
  Q *= 2.5;
  
  for (int i = 0; i <= n; i++)
    if (abs(Q(i) - 2.5*P(i)) > threshold)
      {
	cout << "operateur *= incorrect" << endl;
	abort();
      }

  // on teste l'operateur /=
  Q /= 2.5;
  
  for (int i = 0; i <= n; i++)
    if (abs(Q(i) - P(i)) > threshold)
      {
	cout << "operateur /= incorrect" << endl;
	abort();
      }
  
  // on teste l'operateur +
  for (int i = 0; i <= n; i++)
    P(i) = double(rand())/RAND_MAX;
  
  int m = 7;
  PolynomialTest R(m);
  for (int i = 0; i <= m; i++)
    R(i) = double(rand())/RAND_MAX;
  
  Q = P + R;
  
  for (int i = 0; i < x.GetM(); i++)
    {
      double val = Q.Evaluate(x(i));
      double res = P.Evaluate(x(i)) + R.Evaluate(x(i));
      if (abs(val -res) > threshold)
	{
	  cout << "operateur + incorrect" << endl;
	  abort();
	}
    }

  Q = R + P;
  
  for (int i = 0; i < x.GetM(); i++)
    {
      double val = Q.Evaluate(x(i));
      double res = P.Evaluate(x(i)) + R.Evaluate(x(i));
      if (abs(val -res) > threshold)
	{
	  cout << "operateur + incorrect" << endl;
	  abort();
	}
    }
  
  Q = P + 1.2;
  if (abs(Q(0)-P(0)-1.2) > threshold)
    {
      cout << "operateur + incorrect" << endl;
      abort();
    }
  
  for (int i = 1; i <= n; i++)
    if (Q(i) != P(i))
      {
	cout << "operateur + incorrect" << endl;
	abort();
      }

  Q = 1.2 + P;
  if (abs(Q(0)-P(0)-1.2) > threshold)
    {
      cout << "operateur + incorrect" << endl;
      abort();
    }
  
  for (int i = 1; i <= n; i++)
    if (Q(i) != P(i))
      {
	cout << "operateur + incorrect" << endl;
	abort();
      }

  // on teste l'operateur -
  Q = P - R;
  
  for (int i = 0; i < x.GetM(); i++)
    {
      double val = Q.Evaluate(x(i));
      double res = P.Evaluate(x(i)) - R.Evaluate(x(i));
      if (abs(val -res) > threshold)
	{
	  cout << "operateur - incorrect" << endl;
	  abort();
	}
    }

  Q = R - P;
  
  for (int i = 0; i < x.GetM(); i++)
    {
      double val = Q.Evaluate(x(i));
      double res = R.Evaluate(x(i)) - P.Evaluate(x(i));
      if (abs(val -res) > threshold)
	{
	  cout << "operateur - incorrect" << endl;
	  abort();
	}
    }
  
  Q = P - 1.2;
  if (abs(Q(0)-P(0)+1.2) > threshold)
    {
      cout << "operateur - incorrect" << endl;
      abort();
    }
  
  for (int i = 1; i <= n; i++)
    if (Q(i) != P(i))
      {
	cout << "operateur - incorrect" << endl;
	abort();
      }

  Q = 1.2 - P;
  if (abs(Q(0)+P(0)-1.2) > threshold)
    {
      cout << "operateur - incorrect" << endl;
      abort();
    }
  
  for (int i = 1; i <= n; i++)
    if (abs(Q(i)+P(i)) > threshold)
      {
	cout << "operateur - incorrect" << endl;
	abort();
      }

  Q = -P;
  for (int i = 0; i <= n; i++)
    if (abs(Q(i)+P(i)) > threshold)
      {
	cout << "operateur - incorrect" << endl;
	abort();
      }
  
  // on teste l'operateur *
  Q = P * R;
  
  for (int i = 0; i < x.GetM(); i++)
    {
      double val = Q.Evaluate(x(i));
      double res = P.Evaluate(x(i)) * R.Evaluate(x(i));
      if (abs(val -res) > threshold)
	{
	  cout << "operateur * incorrect" << endl;
	  abort();
	}
    }
  
  Q = P * 1.2;
  for (int i = 0; i <= n; i++)
    if (abs(Q(i) - 1.2*P(i)) > threshold)
      {
	cout << "operateur * incorrect" << endl;
	abort();
      }

  Q = 1.2 * P;
  for (int i = 0; i <= n; i++)
    if (abs(Q(i) - 1.2*P(i)) > threshold)
      {
	cout << "operateur * incorrect" << endl;
	abort();
      }
  
  // on teste la fonction Pow
  Q = Pow(P, 3);
  
  for (int i = 0; i < x.GetM(); i++)
    {
      double val = Q.Evaluate(x(i));
      double res = pow(P.Evaluate(x(i)), 3.0);
      if (abs(val -res) > threshold)
	{
	  cout << "Pow incorrect" << endl;
	  abort();
	}
    }
  
  // on teste la fonction DerivatePolynomial
  DerivatePolynomial(P, Q);
  Real_wp h = 1e-6;
  for (int i = 0; i < x.GetM(); i++)
    {
      double val = Q.Evaluate(x(i));
      double res = (P.Evaluate(x(i)+h) - P.Evaluate(x(i)-h)) / (2.0*h);
      if (abs(val -res) > h)
	{
	  cout << "DerivatePolynomial incorrect" << endl;
	  abort();
	}
    }
  
  // on teste la fonction IntegratePolynomial
  IntegratePolynomial(P, Q);
  for (int i = 0; i < x.GetM(); i++)
    {
      double val = P.Evaluate(x(i));
      double res = (Q.Evaluate(x(i)+h) - Q.Evaluate(x(i)-h)) / (2.0*h);
      if (abs(val -res) > h)
	{
	  cout << "IntegratePolynomial incorrect" << endl;
	  abort();
	}
    }
  
  // on teste la fonction Add
  Q.SetOrder(n);
  for (int i = 0; i <= n; i++)
    Q(i) = double(rand())/RAND_MAX;
  
  P = Q;
  Add(0.45, Plow, P);
  for (int i = 0; i < x.GetM(); i++)
    {
      double val = P.Evaluate(x(i));
      double res = Q.Evaluate(x(i)) + 0.45*Plow.Evaluate(x(i));
      if (abs(val -res) > h)
	{
	  cout << "Add incorrect" << endl;
	  abort();
	}
    }

  P = Q;
  Add(0.52, Phigh, P);
  for (int i = 0; i < x.GetM(); i++)
    {
      double val = P.Evaluate(x(i));
      double res = Q.Evaluate(x(i)) + 0.52*Phigh.Evaluate(x(i));
      if (abs(val -res) > h)
	{
	  cout << "Add incorrect" << endl;
	  abort();
	}
    }

  // on teste la fonction Mlt
  P = Q;
  Mlt(0.28, P);
  for (int i = 0; i < x.GetM(); i++)
    {
      double val = P.Evaluate(x(i));
      double res = 0.28*Q.Evaluate(x(i));
      if (abs(val -res) > h)
	{
	  cout << "Mlt incorrect" << endl;
	  abort();
	}
    }
  
  // on teste la fonction MltAdd
  P = Q;
  MltAdd(0.43, Plow, Phigh, 0.8, P);
  for (int i = 0; i < x.GetM(); i++)
    {
      double val = P.Evaluate(x(i));
      double res = 0.8*Q.Evaluate(x(i)) + 0.43*Plow.Evaluate(x(i))*Phigh.Evaluate(x(i));
      if (abs(val -res) > h)
	{
	  cout << "MltAdd incorrect" << endl;
	  abort();
	}
    }

  // on teste le calcul des polynomes de Legendre
  n = 10;
  Vector<PolynomialTest> Ln(n+1);
  LegendrePolynomial<double> Popt(n+1);
  ComputeLegendrePolynome(0, n+1, Ln);
  Vector<double> PnEval, PnPlus, PnMinus, dPn;
  
  for (int k = 0; k < x.GetM(); k++)
    {
      double y = x(k);
      Popt.EvaluatePn(n+1, y, PnEval);       
      
      double L0 = 1.0, L1 = y;
      for (int i = 0; i <= n; i++)
	{
	  if (isnan(L0) || isnan(Ln(i).Evaluate(y)) || isnan(PnEval(i)) ||
	      (abs(L0 - Ln(i).Evaluate(y)) > threshold) ||
	      (abs(L0 - PnEval(i)) > threshold) )
	    {
	      DISP(i); DISP(y); DISP(L0); DISP(Ln(i).Evaluate(y)); DISP(PnEval(i));
	      cout << "Polynomes de Legendre incorrects" << endl;
	      abort();
	    }            
	  
	  double L2 = ((2.0*i+3)*y*L1 - (i+1)*L0)/double(i+2);
	  L0 = L1;
	  L1 = L2;
	}

      Popt.EvaluatePn(n+1, y+h, PnPlus);
      Popt.EvaluatePn(n+1, y-h, PnMinus);
      Popt.EvaluateDerive(n+1, y, PnEval, dPn);
      for (int i = 0; i <= n; i++)
	{
	  double val = dPn(i), val_ref = (PnPlus(i)-PnMinus(i))/(2.0*h);
	  if (isnan(val) || isnan(val_ref) || abs(val-val_ref) > h)
	    {
	      DISP(i); DISP(val); DISP(val_ref);
	      cout << "EvaluateDerive incorrect" << endl;
	      abort();
	    }
	}
    }

  // on teste le calcul de racines
  Q.SetOrder(n);
  for (int i = 0; i <= n; i++)
    Q(i) = double(rand())/RAND_MAX;

  VectReal_wp Rreal, Rimag;
  SolvePolynomialEquation(Q, Rreal, Rimag);

  for (int i = 0; i < Q.GetOrder(); i++)
    {
      Complex_wp z(Rreal(i), Rimag(i));
      Complex_wp feval = Q.Evaluate(z);
      if (abs(feval) > threshold)
	{
	  cout << "SolvePolynomialEquation" << endl;
	  abort();
	}
    }

  // associated Legendre polynomials are tested
  for (int l = 0; l < 8; l++)
    {
      AssociatedLegendrePolynomial<Real_wp> Pol;
      Pol.Init(l);      
      
      Vector<VectReal_wp> Pnm;
      Real_wp teta = 0.583;
      Pol.EvaluatePnm(l, l, teta, Pnm);

      for (int m = 0; m <= l; m++)
	{
	  R.SetOrder(2);
	  R(0) = Real_wp(-1); R(1) = Real_wp(0); R(2) = Real_wp(1);
	  R = Pow(R, l);
	  for (int j = 0; j < l+m; j++)
	    {
	      DerivatePolynomial(R, Q);
	      R = Q;
	    }

	  // calcul alternatif de Pml
	  Real_wp x = cos(teta);
	  Real_wp val_ref = pow(1.0-x*x, Real_wp(m)/2) / (pow(2.0, Real_wp(l))*fac(l));
	  val_ref *= R.Evaluate(x);
	  Real_wp  coef = sqrt((2*l+1)/(4.0*pi_wp)*fac(l-m)/fac(l+m));
	  if (m%2 == 1)
	    val_ref = -val_ref;

	  if (abs(coef*val_ref - Pnm(l)(m)) > threshold)
	    {
	      DISP(Q);
	      DISP(l); DISP(m); DISP(val_ref); DISP(coef); DISP(coef*val_ref); DISP(Pnm(l)(m));
	      cout << "AssociatedLegendrePolynomial incorrect"  << endl;
	      abort();
	    }
	}
    }
}

template<class PolynomialTest>
void CheckMultivariatePolynomial(PolynomialTest&)
{
  int n = 6;
  MultivariatePolynomial<Real_wp> P(2, n), Q;
  
  if (P.GetM() != (n+1)*(n+2)/2)
    {
      cout << "GetM incorrect" << endl;
      abort();
    }

  if (P.GetDimension() != 2)
    {
      cout << "GetDimension incorrect" << endl;
      abort();
    }
  
  if (P.GetOrder() != n)
    {
      cout << "GetOrder incorrect" << endl;
      abort();
    }

  Q.SetOrder(3, n+2);
  if (Q.GetDimension() != 3)
    {
      cout << "GetDimension/SetOrder incorrect" << endl;
      abort();
    }
  
  if (Q.GetOrder() != n+2)
    {
      cout << "GetOrder/SetOrder incorrect" << endl;
      abort();
    }
  
  int m = Q.GetOrder();
  Array3D<Real_wp> coef(m+1, m+1, m+1);
  coef.Fill(0);
  for (int i = 0; i <= m; i++)
    for (int j = 0; j <= m-i; j++)
      for (int k = 0; k <= m-i-j; k++)
	{
	  coef(i, j, k) = double(rand())/RAND_MAX;
	  Q(i, j, k) = coef(i, j, k);
	}
  
  Vector<double> x;
  Linspace(-1.0, 1.0, 22, x);
  
  VectR3 Points(x.GetM()*x.GetM()*x.GetM());
  int nb = 0;
  for (int ip = 0; ip < x.GetM(); ip++)
    for (int jp = 0; jp < x.GetM(); jp++)
      for (int kp = 0; kp < x.GetM(); kp++)
	{
	  Points(nb).Init(x(ip), x(jp), x(kp));
	  Real_wp val = Q.Evaluate(Points(nb));
	  Real_wp val_ref = 0;
	  for (int i = 0; i <= m; i++)
	    for (int j = 0; j <= m-i; j++)
	      for (int k = 0; k <= m-i-j; k++)
		val_ref += coef(i, j, k)*pow(x(ip), i)*pow(x(jp), j)*pow(x(kp), k);
	  
	  if (isnan(val) || isnan(val_ref) || abs(val-val_ref) > threshold)
	    {
	      DISP(val); DISP(val_ref);
	      cout << "Evaluate incorrect" << endl;
	      abort();
	    }
	  
	  nb++;
	}
  
  Q.ResizeOrder(2, n);
  if (Q.GetDimension() != 2)
    {
      cout << "ResizeOrder incorrect" << endl;
      abort();
    }

  if (Q.GetOrder() != n)
    {
      cout << "ResizeOrder incorrect" << endl;
      abort();
    }
  
  for (int i = 0; i <= n; i++)
    for (int j = 0; j <= n-i; j++)
      {
	if ( (abs(Q(i, j) - coef(i, j, 0)) > threshold) || isnan(Q(i, j)))
	  {
	    DISP(i); DISP(j); DISP(Q(i, j)); DISP(coef(i, j, 0));
	    cout << "ResizeOrder incorrect" << endl;
	    abort();
	  }
      }
  
  Q.ResizeOrder(3, n+2);
  if (Q.GetDimension() != 3)
    {
      cout << "ResizeOrder incorrect" << endl;
      abort();
    }

  if (Q.GetOrder() != n+2)
    {
      cout << "ResizeOrder incorrect" << endl;
      abort();
    }

  for (int i = 0; i <= n; i++)
    for (int j = 0; j <= n-i; j++)
      {
	if ( (abs(Q(i, j, 0) - coef(i, j, 0)) > threshold) || isnan(Q(i, j, 0)))
	  {
	    DISP(i); DISP(j); DISP(Q(i, j, 0)); DISP(coef(i, j, 0));
	    cout << "ResizeOrder incorrect" << endl;
	    abort();
	  }
      }
  
  for (int i = 0; i <= m; i++)
    for (int j = 0; j <= m-i; j++)
      for (int k = 0; k <= m-i-j; k++)
	Q(i, j, k) = coef(i, j, k);
  
  P.Copy(Q);
  for (int np = 0; np < Points.GetM(); np++)
    {
      Real_wp val_ref = Q.Evaluate(Points(np));
      Real_wp val = P.Evaluate(Points(np));
      if (isnan(val) || isnan(val_ref) || abs(val-val_ref) > threshold)
	{
	  DISP(val); DISP(val_ref);
	  cout << "Copy incorrect" << endl;
	  abort();
	}
    }

  P.Fill(0);
  IVect powers(3);
  for (int i = 0; i <= m; i++)
    for (int j = 0; j <= m-i; j++)
      for (int k = 0; k <= m-i-j; k++)
	{
	  powers(0) = i;  powers(1) = j;  powers(2) = k; 
	  int index = P.PowersToIndex(powers);
	  P(index) = coef(i, j, k);
	}
  
  for (int i = 0; i <= m; i++)
    for (int j = 0; j <= m-i; j++)
      for (int k = 0; k <= m-i-j; k++)
	if (abs(P(i, j, k) - coef(i, j, k)) > threshold)
	  {
	    cout << "PowersToIndex incorrect" << endl;
	  }
  
  for (int index = 0; index < P.GetM(); index++)
    {
      P.IndexToPowers(index, powers);
      if (abs(P(index) - coef(powers(0), powers(1), powers(2))) > threshold)
	{
	  cout << "IndexToPowers incorrect" << endl;
	  abort();
	}
    }
  
  P.Clear();
  P = Q;
  for (int np = 0; np < Points.GetM(); np++)
    {
      Real_wp val_ref = Q.Evaluate(Points(np));
      Real_wp val = P.Evaluate(Points(np));
      if (isnan(val) || isnan(val_ref) || abs(val-val_ref) > threshold)
	{
	  DISP(val); DISP(val_ref);
	  cout << "Operator = incorrect" << endl;
	  abort();
	}
    }
  
  // operateur += pour des polynomes
  PolynomialTest Plow(3, m-2), Phigh(3, m+3);
  for (int i = 0; i <= Plow.GetOrder(); i++)
    for (int j = 0; j <= Plow.GetOrder()-i; j++)
      for (int k = 0; k <= Plow.GetOrder()-i-j; k++)
	Plow(i, j, k) = double(rand())/RAND_MAX;

  for (int i = 0; i <= Phigh.GetOrder(); i++)
    for (int j = 0; j <= Phigh.GetOrder()-i; j++)
      for (int k = 0; k <= Phigh.GetOrder()-i-j; k++)
	Phigh(i, j, k) = double(rand())/RAND_MAX;
  
  P = Q;
  P += Plow;
  
  for (int i = 0; i < Points.GetM(); i++)
    {
      double val = P.Evaluate(Points(i));
      double res = Q.Evaluate(Points(i)) + Plow.Evaluate(Points(i));
      if ((abs(val -res) > threshold) || isnan(val) || isnan(res))
	{
	  cout << "operateur += incorrect" << endl;
	  abort();
	}
    }

  P = Q;
  P += Phigh;
  
  for (int i = 0; i < Points.GetM(); i++)
    {
      double val = P.Evaluate(Points(i));
      double res = Q.Evaluate(Points(i)) + Phigh.Evaluate(Points(i));
      if ((abs(val -res) > threshold) || isnan(val) || isnan(res))
	{
	  cout << "operateur += incorrect" << endl;
	  abort();
	}
    }

  // operateur -= pour des polynomes
  P = Q;
  P -= Plow;
  
  for (int i = 0; i < Points.GetM(); i++)
    {
      double val = P.Evaluate(Points(i));
      double res = Q.Evaluate(Points(i)) - Plow.Evaluate(Points(i));
      if ((abs(val -res) > threshold) || isnan(val) || isnan(res))
	{
	  cout << "operateur -= incorrect" << endl;
	  abort();
	}
    }

  P = Q;
  P -= Phigh;
  
  for (int i = 0; i < Points.GetM(); i++)
    {
      double val = P.Evaluate(Points(i));
      double res = Q.Evaluate(Points(i)) - Phigh.Evaluate(Points(i));
      if ((abs(val -res) > threshold) || isnan(val) || isnan(res))
	{
	  cout << "operateur -= incorrect" << endl;
	  abort();
	}
    }
  
  // operateur *= pour des polynomes  
  P = Q;
  P *= Plow;
  
  for (int i = 0; i < Points.GetM(); i++)
    {
      double val = P.Evaluate(Points(i));
      double res = Q.Evaluate(Points(i)) * Plow.Evaluate(Points(i));
      if ((abs(val -res) > 1e3*threshold) || isnan(val) || isnan(res))
	{
	  cout << "operateur *= incorrect" << endl;
	  abort();
	}
    }

  P = Q;
  P *= Phigh;
  
  for (int i = 0; i < Points.GetM(); i++)
    {
      double val = P.Evaluate(Points(i));
      double res = Q.Evaluate(Points(i)) * Phigh.Evaluate(Points(i));
      if ((abs(val -res) > 1e3*threshold) || isnan(val) || isnan(res))
	{
	  cout << "operateur *= incorrect" << endl;
	  abort();
	}
    }
  
  P = Q;
  
  // on teste l'operateur +=
  Q += 0.5;
  if (abs(Q(0, 0, 0)-P(0, 0, 0)-0.5) > threshold)
    {
      cout << "operateur += incorrect" << endl;
      abort();
    }
  
  for (int i = 1; i < Q.GetM(); i++)
    if (Q(i) != P(i))
      {
	cout << "operateur += incorrect" << endl;
	abort();
      }
  
  // on teste l'operateur -=
  Q -= 0.2;
  if (abs(Q(0)-P(0)-0.3) > threshold)
    {
      cout << "operateur -= incorrect" << endl;
      abort();
    }
  
  for (int i = 1; i < Q.GetM(); i++)
    if (Q(i) != P(i))
      {
	cout << "operateur -= incorrect" << endl;
	abort();
      }
  
  // on teste l'operateur *=
  Q(0) = P(0);
  Q *= 2.5;
  
  for (int i = 0; i < Q.GetM(); i++)
    if (abs(Q(i) - 2.5*P(i)) > threshold)
      {
	cout << "operateur *= incorrect" << endl;
	abort();
      }

  // on teste l'operateur +
  P.SetOrder(3, n);
  for (int i = 0; i < P.GetM(); i++)
    P(i) = double(rand())/RAND_MAX;
    
  int mr = 7;
  PolynomialTest R(3, mr);
  for (int i = 0; i < R.GetM(); i++)
    R(i) = double(rand())/RAND_MAX;
  
  Q = P + R;
  
  for (int i = 0; i < Points.GetM(); i++)
    {
      double val = Q.Evaluate(Points(i));
      double res = P.Evaluate(Points(i)) + R.Evaluate(Points(i));
      if (abs(val -res) > threshold)
	{
	  cout << "operateur + incorrect" << endl;
	  abort();
	}
    }

  Q = R + P;
  
  for (int i = 0; i < Points.GetM(); i++)
    {
      double val = Q.Evaluate(Points(i));
      double res = P.Evaluate(Points(i)) + R.Evaluate(Points(i));
      if (abs(val -res) > threshold)
	{
	  cout << "operateur + incorrect" << endl;
	  abort();
	}
    }
  
  Q = P + 1.2;
  if (abs(Q(0)-P(0)-1.2) > threshold)
    {
      cout << "operateur + incorrect" << endl;
      abort();
    }
  
  for (int i = 1; i < Q.GetM(); i++)
    if (Q(i) != P(i))
      {
	cout << "operateur + incorrect" << endl;
	abort();
      }

  Q = 1.2 + P;
  if (abs(Q(0)-P(0)-1.2) > threshold)
    {
      cout << "operateur + incorrect" << endl;
      abort();
    }
  
  for (int i = 1; i < Q.GetM(); i++)
    if (Q(i) != P(i))
      {
	cout << "operateur + incorrect" << endl;
	abort();
      }

  // on teste l'operateur -
  Q = P - R;
  
  for (int i = 0; i < Points.GetM(); i++)
    {
      double val = Q.Evaluate(Points(i));
      double res = P.Evaluate(Points(i)) - R.Evaluate(Points(i));
      if (abs(val -res) > threshold)
	{
	  cout << "operateur - incorrect" << endl;
	  abort();
	}
    }

  Q = R - P;
  
  for (int i = 0; i < Points.GetM(); i++)
    {
      double val = Q.Evaluate(Points(i));
      double res = R.Evaluate(Points(i)) - P.Evaluate(Points(i));
      if (abs(val -res) > threshold)
	{
	  cout << "operateur - incorrect" << endl;
	  abort();
	}
    }
  
  Q = P - 1.2;
  if (abs(Q(0)-P(0)+1.2) > threshold)
    {
      cout << "operateur - incorrect" << endl;
      abort();
    }
  
  for (int i = 1; i < Q.GetM(); i++)
    if (Q(i) != P(i))
      {
	cout << "operateur - incorrect" << endl;
	abort();
      }

  Q = 1.2 - P;
  if (abs(Q(0)+P(0)-1.2) > threshold)
    {
      cout << "operateur - incorrect" << endl;
      abort();
    }
  
  for (int i = 1; i < Q.GetM(); i++)
    if (abs(Q(i)+P(i)) > threshold)
      {
	cout << "operateur - incorrect" << endl;
	abort();
      }

  Q = -P;
  for (int i = 0; i < Q.GetM(); i++)
    if (abs(Q(i)+P(i)) > threshold)
      {
	cout << "operateur - incorrect" << endl;
	abort();
      }
  
  // on teste l'operateur *
  Q = P * R;
  
  for (int i = 0; i < Points.GetM(); i++)
    {
      double val = Q.Evaluate(Points(i));
      double res = P.Evaluate(Points(i)) * R.Evaluate(Points(i));
      if (abs(val -res) > 1e3*threshold)
	{
	  cout << "operateur * incorrect" << endl;
	  abort();
	}
    }
  
  Q = P * 1.2;
  for (int i = 0; i < Q.GetM(); i++)
    if (abs(Q(i) - 1.2*P(i)) > threshold)
      {
	cout << "operateur * incorrect" << endl;
	abort();
      }

  Q = 1.2 * P;
  for (int i = 0; i < Q.GetM(); i++)
    if (abs(Q(i) - 1.2*P(i)) > threshold)
      {
	cout << "operateur * incorrect" << endl;
	abort();
      }
  
  if (P == Q)
    {
      cout << "Operator == incorrect" << endl;
      abort();
    }
  
  Q = P;

  if (P != Q)
    {
      cout << "Operator != incorrect" << endl;
      abort();
    }
  
  // testing function Pow
  P.ResizeOrder(3, n-2);
  
  R = Pow(P, 3);
  for (int i = 0; i < Points.GetM(); i++)
    {
      double val = R.Evaluate(Points(i));
      double res = pow(P.Evaluate(Points(i)), 3.0);
      if ((abs(val -res) > 1e3*threshold) || isnan(val) || isnan(res))
	{
	  cout << "Pow incorrect" << endl;
	  abort();
	}
    }
  
  // testing DerivatePolynomial/IntegratePolynomial
  Real_wp h = 1e-6;
  for (int num = 0; num < 3; num++)
    {
      DerivatePolynomial(Q, P, num);
      R = Derivate(Q, num);
      for (int i = 0; i < Points.GetM(); i++)
	{
	  R3 pt_plus(Points(i)), pt_minus(Points(i));
	  pt_plus(num) += h; pt_minus(num) -= h;
	  double val = P.Evaluate(Points(i));
	  double valr = R.Evaluate(Points(i));
	  double res = (Q.Evaluate(pt_plus) - Q.Evaluate(pt_minus)) / (2.0*h);
	  if ((abs(val -res) > h) || (abs(valr -res) > h) || isnan(val) || isnan(res))
	    {
	      cout << "Derivate incorrect" << endl;
	      abort();
	    }
	}
      
      /*IntegratePolynomial(Q, P, num);
      for (int i = 0; i < Points.GetM(); i++)
	{
	  R3 pt_plus(Points(i)), pt_minus(Points(i));
	  pt_plus(num) += h; pt_minus(num) -= h;
	  double val = Q.Evaluate(Points(i));
	  double res = (P.Evaluate(pt_plus) - P.Evaluate(pt_minus)) / (2.0*h);
	  if ((abs(val -res) > h) || isnan(val) || isnan(res))
	    {
	      cout << "IntegratePolynomial incorrect" << endl;
	      abort();
	    }
	    }      */
    }
  
  // testing curl and div
  TinyVector<MultivariatePolynomial<Real_wp>, 3> Pvec, Qvec;
  Pvec(0).ResizeOrder(3, n); Pvec(1).ResizeOrder(3, n); Pvec(2).ResizeOrder(3, n);
  for (int i = 0; i < Pvec(0).GetM(); i++)
    {
      Pvec(0)(i) = double(rand())/RAND_MAX;
      Pvec(1)(i) = double(rand())/RAND_MAX;
      Pvec(2)(i) = double(rand())/RAND_MAX;
    }
  
  GetCurlPolynomial(Pvec, Qvec);
  GetDivPolynomial(Pvec, P);
  
  // testing Mlt, MltAdd, Add
  P = Q;
  Mlt(1.34, P);
  for (int i = 0; i < Q.GetM(); i++)
    if (isnan(P(i)) || (abs(P(i) - 1.34*Q(i)) > threshold))
      {
	cout << "Mlt incorrect" << endl;
	abort();
      }
  
  P = Q;
  MltAdd(0.38, Plow, Phigh, 0.9, P);
  
  for (int i = 0; i < Points.GetM(); i++)
    {
      double val = P.Evaluate(Points(i));
      double res = 0.9*Q.Evaluate(Points(i)) + 0.38*Plow.Evaluate(Points(i))*Phigh.Evaluate(Points(i));
      if ((abs(val -res) > 1e3*threshold) || isnan(val) || isnan(res))
	{
	  cout << "MltAdd incorrect" << endl;
	  abort();
	}
    }
  
  P = Q;
  Add(0.47, Plow, P);
  for (int i = 0; i < Points.GetM(); i++)
    {
      double val = P.Evaluate(Points(i));
      double res = Q.Evaluate(Points(i)) + 0.47*Plow.Evaluate(Points(i));
      if ((abs(val -res) > threshold) || isnan(val) || isnan(res))
	{
	  cout << "Add incorrect" << endl;
	  abort();
	}
    }

  P = Q;
  Add(0.47, Phigh, P);
  for (int i = 0; i < Points.GetM(); i++)
    {
      double val = P.Evaluate(Points(i));
      double res = Q.Evaluate(Points(i)) + 0.47*Phigh.Evaluate(Points(i));
      if ((abs(val -res) > threshold) || isnan(val) || isnan(res))
	{
	  cout << "Add incorrect" << endl;
	  abort();
	}
    }
  
  // testing IntegrateOnSimplex
  Real_wp value_int = IntegrateOnSimplex(P);
  Real_wp value_ref = 0;
  for (int i = 0; i <= P.GetOrder(); i++)
    for (int j = 0; j <= P.GetOrder()-i; j++)
      for (int k = 0; k <= P.GetOrder()-i-j; k++)
	value_ref += P(i, j, k)*tgamma(i+1)*tgamma(j+1)*tgamma(k+1)/tgamma(i+j+k+4);
  
  DISP(value_int); DISP(value_ref);
  if (isnan(value_int) || isnan(value_ref) || abs(value_int -value_ref) > threshold)
    {
      cout << "IntegrateOnSimplex incorrect" << endl;
      abort();
    }
  
  // testing ScalePolynomVariable
  P = Q;
  Real_wp scale = 0.24;
  ScalePolynomVariable(P, scale);
  for (int i = 0; i < Points.GetM(); i++)
    {
      double val = P.Evaluate(Points(i));
      double res = Q.Evaluate(scale*Points(i));
      if ((abs(val -res) > threshold) || isnan(val) || isnan(res))
	{
	  cout << "ScalePolynomVariable incorrect" << endl;
	  abort();
	}
    }
  
  // testing CompressPolynom
  for (int i = 0; i <= P.GetOrder(); i++)
    for (int j = 0; j <= P.GetOrder()-i; j++)
      for (int k = 0; k <= P.GetOrder()-i-j; k++)
	{
	  P(i, j, k) = 2.0*Real_wp(rand())/RAND_MAX - 1.0;
	  if (i+j+k > n-2)
	    P(i, j, k) *= h;
	}
  
  //DISP(P);
  R = P;
  CompressPolynom(R, h);
  
  if (R.GetOrder() != n-2)
    {
      cout << "CompressPolynom incorrect" << endl;
      abort();
    }
  
  for (int i = 0; i <= R.GetOrder(); i++)
    for (int j = 0; j <= R.GetOrder()-i; j++)
      for (int k = 0; k <= R.GetOrder()-i-j; k++)
	if (abs(P(i, j, k) - R(i, j, k)) > threshold)
	  {
	    cout << "CompressPolynom incorrect" << endl;
	    abort();
	  }
  
  // testing OppositeX, OppositeY, PermuteXY
  R = P;
  OppositeX(R);
  
  for (int i = 0; i < Points.GetM(); i++)
    {
      double val = P.Evaluate(R3(-Points(i)(0), Points(i)(1), Points(i)(2)));
      double res = R.Evaluate(Points(i));
      if ((abs(val -res) > threshold) || isnan(val) || isnan(res))
	{
	  cout << "OppositeX incorrect" << endl;
	  abort();
	}
    }

  R = P;
  OppositeY(R);
  
  for (int i = 0; i < Points.GetM(); i++)
    {
      double val = P.Evaluate(R3(Points(i)(0), -Points(i)(1), Points(i)(2)));
      double res = R.Evaluate(Points(i));
      if ((abs(val -res) > threshold) || isnan(val) || isnan(res))
	{
	  cout << "OppositeX incorrect" << endl;
	  abort();
	}
    }

  R = P;
  PermuteXY(R);
  
  for (int i = 0; i < Points.GetM(); i++)
    {
      double val = P.Evaluate(R3(Points(i)(1), Points(i)(0), Points(i)(2)));
      double res = R.Evaluate(Points(i));
      if ((abs(val -res) > threshold) || isnan(val) || isnan(res))
	{
	  cout << "OppositeX incorrect" << endl;
	  abort();
	}
    }
  
  // calling other functions
  MultivariatePolynomial<Real_wp> S, T, U, V, W;
  GenerateSymPol(P, Q, R, S);
  GenerateSym4_Vertices(P, Q, R, S);
  GenerateSym4_Edges(P, Q, R, S);
  GenerateSym8(P, Q, R, S, T, U, V, W);

}

int main(int argc, char**argv)
{
  InitMontjoie(argc, argv);
  
  UnivariatePolynomial<Real_wp> P;
  
  CheckPolynomial(P);

  MultivariatePolynomial<Real_wp> Q;
  
  CheckMultivariatePolynomial(Q);
  
  cout << "All tests passed successfully" << endl;
  
  return FinalizeMontjoie();
}
