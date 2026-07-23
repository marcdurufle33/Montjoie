//#include "Solver/MontjoieSolver.hxx"
#define MONTJOIE_WITH_THREE_DIM

#include "Mesh/MontjoieMesh.hxx"

using namespace Montjoie;

Real_wp GetRand()
{
  return Real_wp(rand())/RAND_MAX;
}


int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);
  srand(time(NULL));
 
  CylinderParam<double> fct;
  
  R3 ptC;
  R3 vec_u(0, 0, 1);
  R3 vec_v(1, 0, 0), vec_w(0, 1, 0);
  Real_wp radius = 2.8;

  fct.Points.Reallocate(fct.nb_points);
  fct.Normales.Reallocate(fct.nb_points);
  for (int i = 0; i < fct.nb_points; i++)
    {
      Real_wp teta = GetRand()*2.0*pi_wp;
      Real_wp alpha = 10.0*GetRand();
      Real_wp cos_teta = cos(teta), sin_teta = sin(teta);
      
      fct.Points(i) = ptC + alpha*vec_u + radius*(cos_teta*vec_v + sin_teta*vec_w);
      fct.Normales(i) = cos_teta*vec_v + sin_teta*vec_w;
    }

  VectReal_wp xsol(fct.GetN());
  
#ifdef SELDON_WITH_MKL  
  xsol.Zero();
  SolveLeastSquaresMkl(fct, xsol);
  cout << "Solution given by MKL = " << xsol << endl;
#endif

  xsol.Zero();
  SolveLeastSquaresLvm(fct, xsol);
  cout << "Solution given by Levenberg-Marquardt = " << xsol << endl;


  ConicParam<double> fctC;

  Real_wp angle = 0.5*GetRand()*pi_wp;
  cout << "Testing with a cone and angle = " << angle << " Tan = " << tan(angle) << endl;
  
  fctC.Points.Reallocate(fctC.nb_points);
  fctC.Normales.Reallocate(fctC.nb_points);
  for (int i = 0; i < fctC.nb_points; i++)
    {
      Real_wp teta = GetRand()*2.0*pi_wp;
      Real_wp alpha = 10.0*GetRand();
      Real_wp cos_teta = cos(teta), sin_teta = sin(teta);
      Real_wp radius = tan(angle)*alpha;

      fctC.Points(i) = ptC + alpha*vec_u + radius*(cos_teta*vec_v + sin_teta*vec_w);
      R3 vec_plane = cos_teta*vec_v + sin_teta*vec_w;
      fctC.Normales(i) = vec_plane*cos(angle) - vec_u*sin(angle); 
    }
  
  xsol.Reallocate(fctC.GetN());
  xsol(0) = ptC(0); xsol(2) = ptC(1); xsol(2) = ptC(2);
  xsol(3) = vec_u(0); xsol(4) = vec_u(1); xsol(5) = vec_u(2);
  xsol(6) = angle;

  VectReal_wp fvec(fctC.GetM());
  fctC.EvaluateF(xsol, fvec);

#ifdef SELDON_WITH_MKL
  xsol.Zero();
  SolveLeastSquaresMkl(fctC, xsol);

  fctC.EvaluateF(xsol, fvec);

  cout << "Solution given by MKL = " << xsol << endl;
  cout << "Residual = " << Norm2(fvec) << endl;
#endif

  xsol.Zero();
  SolveLeastSquaresLvm(fctC, xsol);
  cout << "Solution given by Levenberg-Marquardt = " << xsol << endl;

  fctC.EvaluateF(xsol, fvec);
  cout << "Residual = " << Norm2(fvec) << endl;
 
  return FinalizeMontjoie();
}

