#include "Solver/MontjoieSolver.hxx"

#include "Elliptic/Maxwell/PhysicalConstant.hxx"
#include "Elliptic/Maxwell/PhysicalConstant.cxx"

using namespace Montjoie;

class NonLinearProblem
{
public :
  Real_wp gammaNL, eps_inf, invEps;
  VectReal_wp source, diag;
  
  void ComputeRightHandSide(const VectReal_wp& En)
  {
    int N = En.GetM();
    source.Reallocate(N);
    diag.Reallocate(N);
    diag.Fill(invEps);
    for (int i = 0; i < N; i++)
      source(i) = En(i)*(eps_inf + En(i)*En(i));        
  }
  
  void ComputeScheme(const VectReal_wp& En, VectReal_wp& scheme)
  {
    for (int i = 0; i < En.GetM(); i++)
      scheme(i) = En(i)*(eps_inf + En(i)*En(i)) - source(i);        
  }
  
  Real_wp GetNorm2Vector(const VectReal_wp& scheme)
  {
    return Norm2(scheme);
  }
  
  void ComputeAndFactoriseDiff(const VectReal_wp& En, VectReal_wp& scale)
  {
    for (int i = 0; i < En.GetM(); i++)
      diag(i) = 1.0/(eps_inf + 3.0*En(i)*En(i));    
  }
  
  void SolveDifferential(const VectReal_wp& En, VectReal_wp& scheme)
  {
    for (int i = 0; i < En.GetM(); i++)
      scheme(i) = En(i)*invEps;
      //scheme(i) = En(i)*diag(i);
  }
  
};

int main(int argc, char** argv)
{
  /* cout << numeric_limits<clock_t>::max() << endl;
  cout << numeric_limits<clock_t>::min() << endl;
  
  clock_t t0 = 2147483646;
  clock_t t1 = -2147483645;
  double diff = double(t1) - double(2.0*numeric_limits<clock_t>::min()) - double(t0);
  cout << "diff = " << diff << endl;
  cout << "arrondi = " << toInteger(round(diff)) << endl;
  int other_choice = 1;
  other_choice += numeric_limits<clock_t>::max() - t0;
  other_choice += (t1 - numeric_limits<clock_t>::min());
  cout << "other solution = " << other_choice << endl;
  
  t0 = -2147483645;
  t1 = -2147483640;
  cout << "diff = " << t1-t0 << endl;
  */

  if (argc < 3)
    {
      cout << "Enter more arguments " << endl;
      abort();
    }
  
  int nb_iter = atoi(argv[1]);
  int nb_points = atoi(argv[2]);
  Real_wp gammaNL = 1e-33;
  Real_wp E0_adim = sqrt(PhysicalConstant::epsilon0_permittivity/gammaNL);
  gammaNL = 1.0;
  Real_wp eps_inf = 1.5;

  Vector<double> En(nb_points), En_next(nb_points);
  Real_wp invEps = 1.0/eps_inf;
  
  NewtonSolver<Real_wp> newton_solver;
  NonLinearProblem var;
  newton_solver.threshold_norme_scheme = 0.01;
  newton_solver.force_reevaluation_jacobian = false;
  newton_solver.print_level = 0;
  var.gammaNL = gammaNL; var.eps_inf = eps_inf; var.invEps = invEps;
  newton_solver.Init(var, En);
  
  for (int k = 0; k < nb_iter; k++)
    {
      En.FillRand();
      Mlt(10.0/E0_adim, En);
      if (false)
	{
	  int n; Real_wp err, x, y, yprime, xsol;
	  for (int i = 0; i < nb_points; i++)
	    {
	      xsol = En(i);
	      Real_wp beta = xsol*(eps_inf + xsol*xsol);
	      
	      // testing Newton algorithm
	      // on the equation alpha x^3 + eps_inf x = beta
	      n = 0; err = 1.0;
	      
	      x = beta*invEps;
	      //y = x*(eps_inf + gammaNL*x*x) - beta;
	      y = x*(eps_inf + x*x) - beta;
	      err = abs(y);
	      while ((n < 10) && (err > 1e-14))
		{
		  yprime = eps_inf + 3.0*x*x;		  
		  x = x - y / yprime;

		  //yprime = eps_inf;
		  //x = x - y*invEps;
		  
		  y = x*(eps_inf + x*x) - beta;
		  err = abs(y);
		  //cout << "Erreur = " << err << endl;
		  n++;
		}
	      
	      //DISP(x);
	      En_next(i) = x;
	    }
	}
      else if (true)
	{
	  Real_wp err, x, y, yprime, xsol, beta;
	  for (int i = 0; i < nb_points; i++)
	    {
	      xsol = En(i);
	      beta = xsol*(eps_inf + xsol*xsol);
	      
	      // testing Newton algorithm
	      // on the equation alpha x^3 + eps_inf x = beta
	      x = beta*invEps;
	      y = x*(eps_inf + x*x) - beta;
	      
	      yprime = eps_inf + 3.0*x*x;
	      x = x - y/yprime;
	      
	      y = x*(eps_inf + x*x) - beta;
	      yprime = eps_inf + 3.0*x*x;
	      x = x - y/yprime;

	      //y = x*(eps_inf + x*x) - beta;
	      //yprime = eps_inf + 3.0*x*x;
	      //x = x - y/yprime;
 
	      //y = x*(eps_inf + x*x) - beta;
	      //if (abs(y) > 1e-14)
	      //{
	      //  cout << "on passe la " << endl;
	      //  yprime = eps_inf + 3.0*x*x;
	      //  x = x - y/yprime;
	      //}
	      
	      //y = x*(eps_inf + x*x) - beta;
	      //err  = abs(y);
	      //cout << "Erreur = " << err << endl;
	      
	      En_next(i) = x;
	    } 
	}
      else if (false)
	{
	  // trying Cardan's formula
	  Real_wp err, x, y, xsol, q_prime, p_prime = eps_inf/3.0, delta_prime;
	  Real_wp one_third = 1.0/3, sqrt_delta, u3, v3;
	  for (int i = 0; i < nb_points; i++)
	    {
	      xsol = En(i);
	      q_prime = -0.5*xsol*(eps_inf + xsol*xsol);
	      // solving x^3 + p x + q = 0
	      delta_prime = -(p_prime*p_prime*p_prime + q_prime*q_prime);
	      sqrt_delta = sqrt(-delta_prime);
	      u3 = -q_prime + sqrt_delta;
	      v3 = -q_prime - sqrt_delta;
	      if (u3 < 0)
		x = -pow(-u3, one_third);
	      else
		x = pow(u3, one_third);

	      if (v3 < 0)
		x += -pow(-v3, one_third);
	      else
		x += pow(v3, one_third);
	      
	      //y = x*(eps_inf + x*x) + 2.0*q_prime;
	      //err  = abs(y);
	      //cout << "Erreur = " << err << endl;
	      En_next(i) = x;
	    }
	}
      else
	{
	  var.ComputeRightHandSide(En);
	  
	  En_next.Fill(0);
	  newton_solver.Solve(var, En_next);
	}
    }
  
  return 0;
}
