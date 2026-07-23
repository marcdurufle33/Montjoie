#include "Output/MontjoieOutput.hxx"

using namespace Montjoie;

template<class T>
void GetRandNumber(T& x)
{
  x = T(rand())/RAND_MAX;
}

template<class T>
void GenerateRandomVector(Vector<T>& x, int n)
{
  x.Reallocate(n);
  for (int i = 0; i < n; i++)
    GetRandNumber(x(i));
}


int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);

  int N = 200;
  VectReal_wp xdiv, ydiv;
  Linspace(0.0, pi_wp, N, xdiv);
  
  ydiv.Reallocate(N);
  for (int i = 0; i < N; i++)
    ydiv(i) = sin(xdiv(i));
  
  int nb_knots = 10;

  VectReal_wp knots;
  Linspace(0.0, pi_wp, nb_knots, knots);
  /* GenerateRandomVector(knots, nb_knots);
  Sort(knots); Mlt(pi_wp, knots);
  knots(0) = 0; knots(nb_knots-1) = pi_wp; */

  VectReal_wp knots_raff;
  Linspace(0.0, pi_wp, 2*(nb_knots-1)+1, knots_raff);
  /* int nb_knots_raff = 2*(nb_knots-1) + 1;
  GenerateRandomVector(knots_raff, nb_knots_raff);
  Sort(knots_raff); Mlt(pi_wp, knots_raff);
  knots_raff(0) = 0; knots_raff(nb_knots_raff-1) = pi_wp;
  
  DISP(knots);
  DISP(knots_raff); */
  
  VectReal_wp xinterp, yinterp, yexact, yinterp_raff, dy_interp, dy_interp_raff, dy_exact;
  int Ninterp = 1000;
  Linspace(0.0, pi_wp, Ninterp, xinterp);

  yexact.Reallocate(Ninterp);
  dy_exact.Reallocate(Ninterp);
  for (int i = 0; i < Ninterp; i++)
    {
      yexact(i) = sin(xinterp(i));
      dy_exact(i) = cos(xinterp(i));
    }

  // we check Init + Evaluate
  for (int order = 1; order <= 6; order++)
    {
      cout << endl << endl << "order = " << order << endl;
      BSplineInterpolation<Real_wp> spline;
      spline.Init(xdiv, ydiv, order, knots);
      
      yinterp.Reallocate(Ninterp);
      for (int i = 0; i < Ninterp; i++)
	yinterp(i) = spline.Evaluate(xinterp(i));
      
      Real_wp err1 = Norm2(yinterp - yexact) / Norm2(yexact);
      cout << "Erreur N = " << err1 << endl;

      BSplineInterpolation<Real_wp> spline_raff;
      spline_raff.Init(xdiv, ydiv, order, knots_raff);
  
      yinterp_raff.Reallocate(Ninterp);
      for (int i = 0; i < Ninterp; i++)
	yinterp_raff(i) = spline_raff.Evaluate(xinterp(i));
      
      Real_wp err2 = Norm2(yinterp_raff - yexact) / Norm2(yexact);
      cout << "Erreur 2N = " << err2 << endl;

      Real_wp order_measured = (log(err1) - log(err2))/log(Real_wp(2));
      cout << "Order of accuracy = " << order_measured << endl;
       if (abs(order_measured - (order+1)) > 0.1)
	{
	  cout << "Init or Evaluate incorrect " << endl;
	  abort();
	}
    }

  // we check EvaluateDerivative
  Real_wp alpha(2.3);
  for (int order = 1; order <= 6; order++)
    {
      cout << endl << endl << "order = " << order << endl;
      BSplineInterpolation<Real_wp> spline;
      spline.Init(xdiv, ydiv, order, knots);
      
      dy_interp.Reallocate(Ninterp);
      yinterp.Reallocate(Ninterp);
      for (int i = 0; i < Ninterp; i++)
	spline.EvaluateDerivative(xinterp(i), yinterp(i), dy_interp(i));
      
      dy_interp.WriteText("dy.dat");
      Real_wp err1 = Norm2(yinterp - yexact) / Norm2(yexact);
      cout << "Erreur N = " << err1 << endl;

      Real_wp errd1 = Norm2(dy_interp - dy_exact) / Norm2(dy_exact);
      cout << "Erreur Derivee N = " << errd1 << endl;

      BSplineInterpolation<Real_wp> spline_raff;
      spline_raff.Init(xdiv, ydiv, order, knots_raff);
  
      yinterp_raff.Reallocate(Ninterp);
      dy_interp_raff.Reallocate(Ninterp);
      for (int i = 0; i < Ninterp; i++)
	spline_raff.EvaluateDerivative(xinterp(i), yinterp_raff(i), dy_interp_raff(i));
      
      Real_wp err2 = Norm2(yinterp_raff - yexact) / Norm2(yexact);
      cout << "Erreur 2N = " << err2 << endl;

      Real_wp errd2 = Norm2(dy_interp_raff - dy_exact) / Norm2(dy_exact);
      cout << "Erreur Derivee N = " << errd1 << endl;

      Real_wp order_measured = (log(err1) - log(err2))/log(Real_wp(2));
      cout << "Order of accuracy = " << order_measured << endl;

      if (abs(order_measured - (order+1)) > 0.1)
	{
	  cout << "Init or EvaluateDerivative incorrect " << endl;
	  abort();
	}
 
      order_measured = (log(errd1) - log(errd2))/log(Real_wp(2));
      cout << "Order of accuracy for y' = " << order_measured << endl;

      if (abs(order_measured - order) > 0.2)
	{
	  cout << "Init or EvaluateDerivative incorrect " << endl;
	  abort();
	}
      
      // we check Mlt
      spline.Mlt(alpha);

      for (int i = 0; i < Ninterp; i++)
	{
	  spline.EvaluateDerivative(xinterp(i), yinterp_raff(i), dy_interp_raff(i));
	  if (abs(alpha*yinterp(i) - yinterp_raff(i)) > 1e-12)
	    {
	      cout << "Mlt incorrect" << endl;
	      abort();
	    }

	  if (abs(alpha*dy_interp(i) - dy_interp_raff(i)) > 1e-12)
	    {
	      cout << "Mlt incorrect" << endl;
	      abort();
	    }
	}

      // we check Clear
      spline.Clear();
    }  

  cout << "All tests passed successfully " << endl;
  
  return FinalizeMontjoie();
}
