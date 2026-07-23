#include "Output/MontjoieOutput.hxx"

using namespace Montjoie;

Real_wp GetRand()
{
  return Real_wp(rand())/RAND_MAX;
}

int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);

  SplineInterpolation<Real_wp> spline;

  int nb_points = 10;
  VectReal_wp x(nb_points), y(nb_points);
  x(0) = 0; x(nb_points-1) = 2.0*pi_wp;
  y(0) = 0; y(nb_points-1) = 0;
  for (int i = 1; i < nb_points-1; i++)
    {
      x(i) = 2.0*pi_wp*GetRand();
      y(i) = sin(x(i));
    }
  
  Sort(nb_points, x, y);
  spline.Init(x, y);

  x.Write("x.dat");
  y.Write("y.dat");

  // checking that first and second-order derivatives are continuous
  Real_wp h = 1e-14;
  for (int i = 1; i < nb_points-1; i++)
    {
      Real_wp y_minus, yp_minus, ypp_minus;
      spline.EvaluateDerivative(x(i)-h, y_minus, yp_minus, ypp_minus);

      Real_wp y_plus, yp_plus, ypp_plus;
      spline.EvaluateDerivative(x(i)+h, y_plus, yp_plus, ypp_plus);
      
      if ((abs(y_minus - y_plus) > 10.0*h) || (abs(yp_minus - yp_plus) > 10.0*h)
	  || (abs(ypp_minus - ypp_plus) > 10.0*h) || isnan(y_minus) || isnan(y_plus)
	  || isnan(yp_minus) || isnan(yp_plus) || isnan(ypp_minus) || isnan(ypp_plus))
	{
	  cout << "Init incorrect" << endl;
	  DISP(i); DISP(y_minus); DISP(y_plus);
	  DISP(yp_minus); DISP(yp_plus);
	  DISP(ypp_minus); DISP(ypp_plus);
	  abort();
	}
    }
  
  // checking that we have a fourth-order interpolation
  int nb_points_interp = 1000;
  VectReal_wp x_interp, y_interp(nb_points_interp);
  Linspace(0.0, 2.0*pi_wp, nb_points_interp, x_interp);
  
  Real_wp err1(0), norme_sol(0);
  for (int i = 0; i < nb_points_interp; i++)
    {
      y_interp(i) = spline.Evaluate(x_interp(i));
      norme_sol += abs(sin(x_interp(i)));
      err1 += abs(y_interp(i) - sin(x_interp(i)));
      if (isnan(y_interp(i)))
	{
	  cout << "Evaluate incorrect" << endl;
	  abort();
	}
    }
 
  err1 /= norme_sol;

  x_interp.Write("x_i.dat");
  y_interp.Write("y_i.dat");
  
  nb_points = 20;
  x.Reallocate(nb_points); y.Reallocate(nb_points);
  x(0) = 0; x(nb_points-1) = 2.0*pi_wp;
  y(0) = 0; y(nb_points-1) = 0;
  for (int i = 1; i < nb_points-1; i++)
    {
      x(i) = 2.0*pi_wp*GetRand();
      y(i) = sin(x(i));
    }
  
  Sort(nb_points, x, y);
  spline.Init(x, y);

  x.Write("x2.dat");
  y.Write("y2.dat");
  
  Real_wp err2(0); norme_sol = Real_wp(0);
  for (int i = 0; i < nb_points_interp; i++)
    {
      y_interp(i) = spline.Evaluate(x_interp(i));
      norme_sol += abs(sin(x_interp(i)));
      err2 += abs(y_interp(i) - sin(x_interp(i)));
      if (isnan(y_interp(i)))
	{
	  cout << "Evaluate incorrect" << endl;
	  abort();
	}
    }
  
  err2 /= norme_sol;

  y_interp.Write("y_i2.dat");
  
  DISP(err1); DISP(err2); DISP(err1/err2);
  if (isnan(err1) || isnan(err1) || (err2 == Real_wp(0)) || abs(err1/err2) < 16)
    {
      cout << "Evaluate incorrect (not fourth order)" << endl;
      abort();
    }
  
  // checking EvaluateDerivative
  h = 1e-5;
  for (int i = 1; i < nb_points_interp-1; i++)
    {
      Real_wp xp = x_interp(i) + h, xm = x_interp(i) -h;
      Real_wp fp = spline.Evaluate(xp);
      Real_wp fm = spline.Evaluate(xm);
      Real_wp f = spline.Evaluate(x_interp(i));
      Real_wp df = (fp - fm)/(2.0*h);
      Real_wp d2f = (fp - 2.0*f + fm)/(h*h);
      Real_wp yi, dyi, d2yi;
      spline.EvaluateDerivative(x_interp(i), yi, dyi, d2yi);
      if (isnan(yi) || isnan(dyi) || isnan(d2yi) ||
	  abs(yi - f) > h*h || abs(dyi - df) > h*h || abs(d2yi-d2f) > h)
	{
	  cout << "EvaluateDerivative incorrect" << endl;
	  DISP(yi); DISP(f); DISP(dyi); DISP(df); DISP(d2yi); DISP(d2f);
	  DISP(i); DISP(x_interp(i)); DISP(yi-f); DISP(dyi-df); DISP(d2yi-d2f);
	  abort();
	}
    }

  // checking that GetSecondDerivatives is correct
  const VectReal_wp& ypp = spline.GetSecondDerivatives();
  for (int i = 0; i < nb_points; i++)
    {
      Real_wp yi, dyi, d2yi;
      spline.EvaluateDerivative(x(i), yi, dyi, d2yi);
      if (isnan(d2yi) || isnan(ypp(i)) || abs(d2yi - ypp(i)) > 1e-12)
	{
	  cout << "GetSecondDerivatives incorrect" << endl;
	  abort();
	}
    }

  // checking that FindInterval is correct
  for (int i = 0; i < nb_points_interp; i++)
    {
      Real_wp xi = GetRand()*2.0*pi_wp;
      int a, b;
      spline.FindInterval(xi, a, b);
      if ((xi < x(a) -1e-12) || (xi > x(b) + 1e-12))
	{
	  cout << "FindInterval incorrect" << endl;
	  DISP(a); DISP(b); DISP(x(a)); DISP(xi);
	  abort();
	}
    }

  // checking EvaluateFctBasis
  SplineInterpolation<Real_wp> spline2;
  for (int i = 0; i < nb_points; i++)
    y(i) = 2.0*GetRand() - 1.0;
  
  spline2.Init(x, y);

  VectReal_wp vec_phi(nb_points);
  for (int i = 0; i < nb_points_interp; i++)
    {
      spline.EvaluateFctBasis(x_interp(i), vec_phi);
      Real_wp fphi = 0;
      for (int j = 0; j < vec_phi.GetM(); j++)
	fphi += vec_phi(j)*y(j);
      
      Real_wp f = spline2.Evaluate(x_interp(i));
      if (isnan(fphi) || isnan(f) || (abs(f-fphi) > 1e-12))
	{
	  cout << "EvaluateFctBasis incorrect" << endl;
	  abort();
	}
    }

  // checking Mlt
  Real_wp alpha(2.3);
  spline.Init(x, y);
  spline.Mlt(alpha);
  for (int i = 0; i < nb_points_interp; i++)
    {
      Real_wp f = spline.Evaluate(x_interp(i));
      Real_wp fref = spline2.Evaluate(x_interp(i));
      
      if (isnan(f) || isnan(fref) || (abs(f-alpha*fref) > 1e-12))
	{
	  cout << "Mlt incorrect" << endl;
	  abort();
	}
    }

  spline.Clear();
 
  return FinalizeMontjoie();
}
