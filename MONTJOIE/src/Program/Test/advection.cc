#define MONTJOIE_WITH_ONE_DIM

#include "Hyperbolic/WaterWaves/MontjoieKdv.hxx"

using namespace Montjoie;

int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);
    
  if (argc != 2)
    {
      cout << "you need to enter a data file : ./kdv.x test.ini" << endl;
      abort();
    }
  
  cerr.precision(15);
  cout.precision(15);
  
  InputVariablesWaterWaves vars;
  ReadInputFile(string(argv[1]), vars);
   
  double dt = vars.dt;
  double Tfinal = vars.Tfinal;
  int nb_iterations = toInteger(ceil(Tfinal/dt));
  
  // computational domain
  int order_time = vars.time_order;
  int type_dg_scheme = TimeSchemeEnum::RUNGE_KUTTA;
  int type_scheme = TimeSchemeEnum::LEAP_FROG;
  
  int Nd = vars.number_grid_points; // number of points where we want to know the solution
  int delta_display = vars.delta_snapshot; 
  
  double freq = vars.frequency, xcenter = vars.center_source;
  int type_init = vars.type_initial_condition;
  int type_bottom = vars.type_bottom;
  
  Real_wp x0 = vars.x0, xN = vars.xN;
  Real_wp xg_0 = vars.xg_0, xg_N = vars.xg_N;
  
  int nb_elt = vars.nb_elt, order = vars.order;  
  
  double mu = vars.mu;
  double epsilon = vars.cte_epsilon*pow(mu, vars.power_epsilon);
  double alpha = vars.cte_alpha*pow(mu, vars.power_alpha);
  double beta = vars.cte_beta*pow(mu, vars.power_beta);
  
  DISP(mu); DISP(epsilon); DISP(alpha); DISP(beta);
  
  // initialisation of the spectral method
  TimeSpectralScheme1D<AdvectionEquationSpectral> sp_scheme;

  sp_scheme.SetInterval(x0, xN, nb_elt, order);
  sp_scheme.InitGrid(xg_0, xg_N, Nd);
  
  sp_scheme.var_eq.InitEpsilon(epsilon, mu, alpha, beta);
  sp_scheme.var_eq.SetInitialCondition(type_init, type_bottom, xcenter, freq);
  
  sp_scheme.SetOrder(order);
  sp_scheme.InitTimeScheme(type_scheme, order_time, dt, Tfinal);
  
  // initialisation of the finite element method
  TimeLdgScheme1D<AdvectionEquationFem> fem_scheme;      
    
  fem_scheme.SetInterval(x0, xN, nb_elt, order);
  fem_scheme.InitGrid(xg_0, xg_N, Nd);
  
  fem_scheme.var_eq.InitEpsilon(epsilon, mu, alpha, beta);
  fem_scheme.var_eq.SetInitialCondition(type_init, type_bottom, xcenter, freq);
  
  fem_scheme.SetOrder(order);
  fem_scheme.InitTimeScheme(type_dg_scheme, order_time, dt, Tfinal);
    
  // initialisation of the finite difference method
  FiniteDifferenceScheme1D<AdvectionEquation> fd_scheme;
  
  fd_scheme.SetInterval(x0, xN, nb_elt, order);
  fd_scheme.InitGrid(xg_0, xg_N, Nd);
  
  fd_scheme.var_eq.InitEpsilon(epsilon, mu, alpha, beta);
  fd_scheme.var_eq.SetInitialCondition(type_init, type_bottom, xcenter, freq);
  
  fd_scheme.SetOrder(order);
  fd_scheme.InitTimeScheme(type_scheme, order_time, dt, Tfinal);
  
  // initialisation of the LDG method
  TimeLdgScheme1D<AdvectionEquationDG> dg_scheme;      
  
  
  dg_scheme.SetInterval(x0, xN, nb_elt, order);
  dg_scheme.InitGrid(xg_0, xg_N, Nd);
  
  dg_scheme.var_eq.InitEpsilon(epsilon, mu, alpha, beta);
  dg_scheme.var_eq.SetInitialCondition(type_init, type_bottom, xcenter, freq);
  
  dg_scheme.SetOrder(order);
  dg_scheme.InitTimeScheme(type_dg_scheme, order_time, dt, Tfinal);
  
  // storing snapshots
  DISP(nb_iterations); DISP(Nd);
  Matrix<double> Un_sp(Nd, nb_iterations/delta_display+1);
  Matrix<double> Un_fem(Nd, nb_iterations/delta_display+1);
  Matrix<double> Un_dg(Nd, nb_iterations/delta_display+1);
  Matrix<double> Un_fd(Nd, nb_iterations/delta_display+1);
  Un_sp.Fill(0);
  Un_fem.Fill(0);
  Un_dg.Fill(0);
  Un_fd.Fill(0);
  Vector<double> Un;
  for (int n = 0; n <= nb_iterations; n++)
    {
      
      if (n%delta_display == 0)
	{
	  if (n%100 == 0)
	    {
	      DISP(n*dt);
	      cout << "Energie dg " << dg_scheme.GetEnergy() << endl;
	    }
	  
	  int nb = n/delta_display;
	  
	  fd_scheme.GetInterpolateUn(0, Un, fd_scheme.GetIterate());
	  for (int i = 0; i < Nd; i++)
	    Un_fd(i, nb) = Un(i);
	  
	  fem_scheme.GetInterpolateUn(0, Un);
	  for (int i = 0; i < Nd; i++)
	    Un_fem(i, nb) = Un(i);
	  
	  dg_scheme.GetInterpolateUn(0, Un);
	  for (int i = 0; i < Nd; i++)
	    Un_dg(i, nb) = Un(i);
	  
	  sp_scheme.GetInterpolateUn(0, Un);	      
	  for (int i = 0; i < Nd; i++)
	    Un_sp(i, nb) = Un(i);	  
	}
      
      fd_scheme.Advance(n*dt, n);
      dg_scheme.Advance(n*dt, n);
      fem_scheme.Advance(n*dt, n);
      sp_scheme.Advance(n*dt, n);
    }
  
  Un_sp.WriteText("UnSpec.dat");    
  Un_fem.WriteText("UnFem.dat");    
  Un_dg.WriteText("UnRef.dat");    
  Un_fd.WriteText("UnDiff.dat");    
  
  return FinalizeMontjoie();
}
