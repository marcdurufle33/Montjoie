#define MONTJOIE_WITH_ONE_DIM
#define MONTJOIE_WITH_NODAL_H1

#include "Hyperbolic/WaterWaves/MontjoieCamassaHolm.hxx"

using namespace Montjoie;

template<class GenericPb>
void RunAll(GenericPb& scheme, InputVariablesWaterWaves& vars)
{
  typedef Montjoie::Real_wp Real_wp;
  
  double dt = vars.dt;
  double Tfinal = vars.Tfinal;
  int nb_iterations = toInteger(ceil(Tfinal/dt));
  
  // computational domain
  int order_time = vars.time_order;
  int type_scheme = vars.type_time_scheme;
  
  int Nd = vars.number_grid_points; // number of points where we want to know the solution
  int delta_display = vars.delta_snapshot; 
  
  double log_eps_min = log10(vars.mu);
  double log_eps_max = log10(vars.mu);
  int nb_points = 1;
  
  if (vars.log_log_curve)
    {
      log_eps_min = vars.log_eps_min;
      log_eps_max = vars.log_eps_max;
      nb_points = vars.nb_points_log_eps;
    }
  
  int nb_points_effective = nb_points;
  Matrix<double> result_error(nb_points_effective, 3);
  result_error.Fill(0);

  double freq = vars.frequency, xcenter = vars.center_source;
  int type_init = vars.type_initial_condition;
  int type_bottom = vars.type_bottom;
  int type_model = vars.type_model;

  Real_wp x0 = vars.x0, xN = vars.xN;
  Real_wp xg_0 = vars.xg_0, xg_N = vars.xg_N;
  
  int nb_elt = vars.nb_elt, order = vars.order;  
  
#ifdef SELDON_WITH_MPI
  int rank_proc; MPI_Comm_rank(MPI_COMM_WORLD, &rank_proc);
  int nb_proc; MPI_Comm_size(MPI_COMM_WORLD, &nb_proc);
  //int nplot = 0;
  nb_points_effective = (nb_points-rank_proc+nb_proc-1)/nb_proc;
  if (rank_proc > 0)
    result_error.Reallocate(nb_points_effective, 3);
      
  for (int num_point = rank_proc; num_point < nb_points; num_point += nb_proc)
#else
    for (int num_point = 0; num_point < nb_points; num_point++)
#endif
      {
	double log_eps = log_eps_min;
	if (nb_points > 1)
	  log_eps = log_eps_min + (log_eps_max - log_eps_min)*double(num_point)/(nb_points-1);
	
	double mu = pow(10, log_eps), epsilon = vars.cte_epsilon*pow(mu, vars.power_epsilon);
        double alpha = vars.cte_alpha*pow(mu, vars.power_alpha);
        double beta = vars.cte_beta*pow(mu, vars.power_beta);
                
        DISP(mu); DISP(epsilon); DISP(alpha); DISP(beta);
	
	cout << "Computing solution for epsilon " << epsilon << endl;
	
	// initialization of the scheme
        scheme.var_eq.type_model = type_model;
	scheme.SetInterval(x0, xN, nb_elt, order);
	scheme.InitGrid(xg_0, xg_N, Nd);
	
	scheme.var_eq.InitEpsilon(epsilon, mu, alpha, beta);
        scheme.var_eq.SetInitialCondition(type_init, type_bottom,
					  xcenter, freq);
        scheme.SetOrder(order);
        scheme.InitTimeScheme(type_scheme, order_time, dt, Tfinal);
	
	//VectReal_wp ProdUn(scheme.GetIterateReal());
	//scheme.EvaluateFunction(0.0, scheme.GetIterateReal(), ProdUn);
	//scheme.GetIterateReal().WriteText("UnInit.dat");
	//ProdUn.WriteText("ProdUnInit.dat");
        
        // storing snapshots
	Matrix<double> Un_history(Nd, nb_iterations/delta_display+1);
	Un_history.Fill(0);
	VectReal_wp Un;
        
	// main loop in time
	for (int n = 0; n <= nb_iterations; n++)
	  {
            if (n%delta_display == 0)
	      {
		DISP(n*dt);
                cout << "Energie du schema " << scheme.GetEnergy() << endl;
		
		int nb = n/delta_display;
                scheme.GetInterpolateUn(0, Un);
                
		for (int i = 0; i < Nd; i++)
		  Un_history(i, nb) = Un(i);
              }
            
            scheme.Advance(n*dt, n);
          }
        
        DISP(scheme.GetNbEvaluationDifferential());
        Un_history.WriteText(vars.file_output_history);
      }
}
  

int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);
    
  if (argc != 2)
    {
      cout << "you need to enter a data file : ./kdv.x test.ini" << endl;
      abort();
    }
  
  cerr.precision(17);
  cout.precision(17);
  
  InputVariablesWaterWaves vars;
  ReadInputFile(string(argv[1]), vars);
  
  if (vars.type_formulation == vars.SPECTRAL)
    {
      cout << "Use of a spectral method to solve Kdv equation" << endl;
      TimeSpectralScheme1D<CamassaHolmEquationSpectral> scheme;
      RunAll(scheme, vars);
    }
  else if (vars.type_formulation == vars.FDTD)
    {
      cout << "Use of a finite-difference method to solve Kdv equation" << endl;
      FiniteDifferenceScheme1D<CamassaHolmEquation> scheme;
      RunAll(scheme, vars);
    }
  else if (vars.type_formulation == vars.LDG)
    {
      cout << "Use of a discontinuous Galerkin method to solve Kdv equation" << endl;
      TimeLdgScheme1D<CamassaHolmEquationDG> scheme;
      RunAll(scheme, vars);
    }
  else if (vars.type_formulation == vars.FEM)
    {
      cout << "Use of a finite element method to solve Kdv equation" << endl;
      TimeLdgScheme1D<CamassaHolmEquationFem> scheme;
      RunAll(scheme, vars);
    }
  
  return FinalizeMontjoie();
}
