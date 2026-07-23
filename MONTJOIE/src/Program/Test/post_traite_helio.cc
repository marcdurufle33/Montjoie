#include "Output/MontjoieOutput.hxx"

using namespace Montjoie;

void RefineVector(const VectReal_wp& x, VectReal_wp& y, int r)
{
  int nb_div = x.GetM() - 1;
  y.Reallocate(nb_div*r + 1);
  for (int i = 0; i < nb_div; i++)
    {
      Real_wp a = x(i), b = x(i+1);
      Real_wp dx = (b-a)/r;
      for (int j = 0; j < r; j++)
	y(i*r + j) = a + j*dx;
    }
  
  y(0) = x(0); y(nb_div*r) = x(nb_div);
}

int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);
  
  if (argc < 4)
    {
      cout << "Entrez un nom de fichier a post-traiter" << endl;
      abort();
    }
  
  string file_name(argv[1]);
  int order = atoi(argv[2]);
  string file_knots(argv[3]);
  bool reverse = true;
  if (argc == 5)
    reverse = false;
  
  Matrix<double> val;
  val.ReadText(file_name);

  int N = val.GetM();
  VectReal_wp rho(N), c(N), r(N), log_rho(N);
  for (int i = 0; i < N; i++) 
    {
      if (reverse)
	{
	  c(i) = val(N-1-i, 1)/100;
	  r(i) = val(N-1-i, 0);
	  rho(i) = val(N-1-i, 2)*1e3;
	}
      else
	{
	  r(i) = val(i, 0);
	  c(i) = val(i, 1);
	  rho(i) = val(i, 2);
	}
      
      log_rho(i) = log(rho(i));
    }
  
  cout << "Last value of c = " << c(N-1) << '\n' << endl;
  
  // technique avec une spline cubique
  /* SplineInterpolation<Real_wp> spline_rho;
  spline_rho.Init(r, rho);
  
  ofstream file_rhop_over_rho("rhop_div_rho.don"), file_velocity("velocity.don");
  ofstream file_radius("radius.don"), file_rho("rho.don");
  file_rhop_over_rho.precision(15); file_velocity.precision(15);
  file_radius.precision(15); file_rho.precision(15);
  Real_wp rho_eval, rhop_eval, rhopp_eval;
  for (int i = 0; i < N; i++) 
    {
      spline_rho.EvaluateDerivative(r(i), rho_eval, rhop_eval, rhopp_eval);
      file_rhop_over_rho << r(i) << " " << rhop_eval / rho_eval << '\n';
      file_velocity << r(i) << " " << c(i) << '\n';
      file_radius << r(i) << '\n';
      file_rho << r(i) << " " << rho(i) << '\n';
    }
  
  file_rho.close();
  file_radius.close();
  file_rhop_over_rho.close();
  file_velocity.close(); */

  // technique avec des B-splines
  DISP(N);

  // on lit les knots dans un fichier
  VectReal_wp knots;
  knots.ReadText(file_knots);
  DISP(knots);
  
  BSplineInterpolation<Real_wp> spline_log_rho;
  spline_log_rho.Init(r, log_rho, order, knots);

  BSplineInterpolation<Real_wp> spline_c0;
  spline_c0.Init(r, c, order, knots);

  VectReal_wp points = spline_log_rho.GetPoints();
  VectReal_wp weights = spline_log_rho.GetWeights();
  VectReal_wp points_c0 = spline_c0.GetPoints();
  VectReal_wp weights_c0 = spline_c0.GetWeights();
 
  //VectReal_wp r_interp;
  //RefineVector(r, r_interp, 4);
  VectReal_wp r_interp;
  Linspace(knots(0), knots(knots.GetM()-1), 20000, r_interp);

  ofstream file_rho_i("rho_interp.don"); file_rho_i.precision(15);
  ofstream file_c0_i("c0_interp.don"); file_c0_i.precision(15);
  for (int i = 0; i < r_interp.GetM(); i++)
    {
      //file_rho_i << r_interp(i) << " " << exp(spline_log_rho.Evaluate(r_interp(i))) << '\n';
      //file_c0_i << r_interp(i) << " " << spline_c0.Evaluate(r_interp(i)) << '\n';

      Real_wp dlog_rho, log_rho, c0, dc0;
      spline_log_rho.EvaluateDerivative(r_interp(i), log_rho, dlog_rho);
      spline_c0.EvaluateDerivative(r_interp(i), c0, dc0);

      file_rho_i << r_interp(i) << " " << exp(log_rho) << " " << exp(log_rho)*dlog_rho << '\n';
      file_c0_i << r_interp(i) << " " << c0 << " " << dc0 << '\n';
    }
  
  file_rho_i.close();
  file_c0_i.close();

  ofstream file_radius(string("radius_r"+to_str(order)+".don").data());
  file_radius.precision(16);
  ofstream file_rho(string("log_rho_r"+to_str(order)+".don").data());
  file_rho.precision(16);
  ofstream file_c0(string("c0_r"+to_str(order)+".don").data());
  file_c0.precision(16);
  for (int i = 0; i < points.GetM(); i++)
    {
      Real_wp poids(0), poids_c0(0);
      if (i < weights.GetM())
	poids = weights(i);

      if (i < weights_c0.GetM())
	poids_c0 = weights_c0(i);
      
      if ((i >= order) && (i < points.GetM()-order))
	file_radius << points(i) << '\n';
      
      file_rho << points(i) << " " << poids << '\n';
      file_c0 << points_c0(i) << " " << poids_c0 << '\n';
    }
  
  file_radius.close();
  file_rho.close();
  file_c0.close();
  
  return FinalizeMontjoie();
}
