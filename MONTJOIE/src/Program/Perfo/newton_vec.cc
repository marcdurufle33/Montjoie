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
    Real_wp En_module;
    for (int i = 0; i < N/2; i++)
      {
	En_module = En(2*i)*En(2*i) + En(2*i+1)*En(2*i+1);
	source(2*i) = En(2*i)*(eps_inf + En_module);
	source(2*i+1) = En(2*i+1)*(eps_inf + En_module);
      }
  }
  
  void ComputeScheme(const VectReal_wp& En, VectReal_wp& scheme)
  {
    Real_wp indice;
    for (int i = 0; i < En.GetM()/2; i++)
      {
	indice = eps_inf + En(2*i)*En(2*i) + En(2*i+1)*En(2*i+1);
	scheme(2*i) = En(2*i)*indice - source(2*i);
	scheme(2*i+1) = En(2*i+1)*indice - source(2*i+1);
      }
  }
  
  Real_wp GetNorm2Vector(const VectReal_wp& scheme)
  {
    return Norm2(scheme);
  }
  
  void ComputeAndFactoriseDiff(const VectReal_wp& En, VectReal_wp& scale)
  {
    Real_wp indice;
    for (int i = 0; i < En.GetM()/2; i++)
      {
	indice = eps_inf + En(2*i)*En(2*i) + En(2*i+1)*En(2*i+1);
	diag(2*i) = 1.0/(eps_inf + indice);
	diag(2*i+1) = 1.0/(eps_inf + indice);
      }    
  }
  
  void SolveDifferential(const VectReal_wp& En, VectReal_wp& scheme)
  {
    for (int i = 0; i < En.GetM(); i++)
      scheme(i) = En(i)*diag(i);
  }
  
};

int main(int argc, char** argv)
{
  cout.precision(15);
  
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

  Vector<double> En(2*nb_points), En_next(2*nb_points);
  Real_wp invEps = 1.0/eps_inf;
  
  NewtonSolver<Real_wp> newton_solver;
  NonLinearProblem var;
  newton_solver.threshold_norme_scheme = 0.01;
  newton_solver.force_reevaluation_jacobian = false;
  newton_solver.print_level = 0;
  newton_solver.nbitermax = 100;
  var.gammaNL = gammaNL; var.eps_inf = eps_inf; var.invEps = invEps;
  newton_solver.Init(var, En);
  
  for (int k = 0; k < nb_iter; k++)
    {
      En.FillRand();
      Mlt(1.0/E0_adim, En);
      if (false)
	{
	  int n; Real_wp err, Ex, Ey,  Ex_sol, Ey_sol;
	  TinyMatrix<Real_wp, Symmetric, 2, 2> mat_DF;
	  TinyVector<Real_wp, 2> vec_F, vec_V;
	  Real_wp beta_x, beta_y, E_module, inv_det;
	  for (int i = 0; i < nb_points; i++)
	    {
	      Ex_sol = En(2*i);
	      Ey_sol = En(2*i+1);
	      E_module = eps_inf + Ex_sol*Ex_sol + Ey_sol*Ey_sol;
	      beta_x = Ex_sol*E_module;
	      beta_y = Ey_sol*E_module;
	      
	      Ex = beta_x*invEps;
	      Ey = beta_y*invEps;
	      
	      E_module = eps_inf + Ex*Ex + Ey*Ey;
	      vec_F(0) = Ex*E_module - beta_x;
	      vec_F(1) = Ey*E_module - beta_y;
	      err = Norm1(vec_F);
	      while ((n < 10) && (err > 1e-14))
		{
		  if (true)
		    {
		      mat_DF(0, 0) = eps_inf + Ey*Ey + 3.0*Ex*Ex;
		      mat_DF(0, 1) = 2.0*Ex*Ey;
		      mat_DF(1, 1) = eps_inf + Ex*Ex + 3.0*Ey*Ey;
		      inv_det = 1.0/Det(mat_DF);
		      
		      vec_V(0) = inv_det*(mat_DF(1, 1)*vec_F(0) - mat_DF(0, 1)*vec_F(1));
		      vec_V(1) = inv_det*(mat_DF(0, 0)*vec_F(1) - mat_DF(0, 1)*vec_F(0));
		    }
		  else
		    {
		      vec_V(0) = invEps*vec_F(0);
		      vec_V(1) = invEps*vec_F(1);
		    }
		  
		  Ex -= vec_V(0);
		  Ey -= vec_V(1);
		  
		  E_module = eps_inf + Ex*Ex + Ey*Ey;
		  vec_F(0) = Ex*E_module - beta_x;
		  vec_F(1) = Ey*E_module - beta_y;
		  err = Norm1(vec_F);
		  //cout << "Residu = " << err << endl;
		  n++;
		}
	      
	      //DISP(Ex_sol); DISP(Ey_sol);
	      //DISP(Ex); DISP(Ey);
	      En_next(2*i) = Ex;
	      En_next(2*i+1) = Ey;
	    }
	}
      else if (true)
	{
	  int n; Real_wp err, Ex, Ey,  Ex_sol, Ey_sol;
	  TinyMatrix<Real_wp, Symmetric, 2, 2> mat_DF;
	  TinyVector<Real_wp, 2> vec_F, vec_V;
	  Real_wp beta_x, beta_y, E_module, inv_det;
	  for (int i = 0; i < nb_points; i++)
	    {
	      Ex_sol = En(2*i);
	      Ey_sol = En(2*i+1);
	      E_module = eps_inf + Ex_sol*Ex_sol + Ey_sol*Ey_sol;
	      beta_x = Ex_sol*E_module;
	      beta_y = Ey_sol*E_module;
	      
	      Ex = beta_x*invEps;
	      Ey = beta_y*invEps;
	      
	      E_module = eps_inf + Ex*Ex + Ey*Ey;
	      vec_F(0) = Ex*E_module - beta_x;
	      vec_F(1) = Ey*E_module - beta_y;
	      
	      mat_DF(0, 0) = eps_inf + Ey*Ey + 3.0*Ex*Ex;
	      mat_DF(0, 1) = 2.0*Ex*Ey;
	      mat_DF(1, 1) = eps_inf + Ex*Ex + 3.0*Ey*Ey;
	      inv_det = 1.0/Det(mat_DF);
	      
	      vec_V(0) = inv_det*(mat_DF(1, 1)*vec_F(0) - mat_DF(0, 1)*vec_F(1));
	      vec_V(1) = inv_det*(mat_DF(0, 0)*vec_F(1) - mat_DF(0, 1)*vec_F(0));
	      
	      Ex -= vec_V(0);
	      Ey -= vec_V(1);
	      
	      E_module = eps_inf + Ex*Ex + Ey*Ey;
	      vec_F(0) = Ex*E_module - beta_x;
	      vec_F(1) = Ey*E_module - beta_y;
	      
	      mat_DF(0, 0) = eps_inf + Ey*Ey + 3.0*Ex*Ex;
	      mat_DF(0, 1) = 2.0*Ex*Ey;
	      mat_DF(1, 1) = eps_inf + Ex*Ex + 3.0*Ey*Ey;
	      inv_det = 1.0/Det(mat_DF);
	      
	      vec_V(0) = inv_det*(mat_DF(1, 1)*vec_F(0) - mat_DF(0, 1)*vec_F(1));
	      vec_V(1) = inv_det*(mat_DF(0, 0)*vec_F(1) - mat_DF(0, 1)*vec_F(0));
	      
	      Ex -= vec_V(0);
	      Ey -= vec_V(1);
	      
	      E_module = eps_inf + Ex*Ex + Ey*Ey;
	      vec_F(0) = Ex*E_module - beta_x;
	      vec_F(1) = Ey*E_module - beta_y;
	      err = Norm1(vec_F);
	      if (err > 1e-14)
		{
		  //cout << "on passe la " << endl;
		  mat_DF(0, 0) = eps_inf + Ey*Ey + 3.0*Ex*Ex;
		  mat_DF(0, 1) = 2.0*Ex*Ey;
		  mat_DF(1, 1) = eps_inf + Ex*Ex + 3.0*Ey*Ey;
		  inv_det = 1.0/Det(mat_DF);
		  
		  vec_V(0) = inv_det*(mat_DF(1, 1)*vec_F(0) - mat_DF(0, 1)*vec_F(1));
		  vec_V(1) = inv_det*(mat_DF(0, 0)*vec_F(1) - mat_DF(0, 1)*vec_F(0));
		  
		  Ex -= vec_V(0);
		  Ey -= vec_V(1);	      
		}
	      //cout << "Residu = " << err << endl;
	      
	      //DISP(Ex_sol); DISP(Ey_sol);
	      //DISP(Ex); DISP(Ey);
	      En_next(2*i) = Ex;
	      En_next(2*i+1) = Ey;
	    }
	}
      else
	{
	  var.ComputeRightHandSide(En);
	  
	  En_next.Fill(0);
	  newton_solver.Solve(var, En_next);
	  //DISP(En); DISP(En_next);
	}
    }
  
  return 0;
}
