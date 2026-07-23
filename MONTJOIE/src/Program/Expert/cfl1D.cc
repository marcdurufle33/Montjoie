#define MONTJOIE_WITH_TWO_DIM

#include "Quadrature/MontjoieQuadrature.hxx"

using namespace Montjoie;

#ifdef MONTJOIE_WITH_MPFR
Real_wp kmax(0,MONTJOIE_DEFAULT_PRECISION), threshold_rational(0,MONTJOIE_DEFAULT_PRECISION);
#else
Real_wp kmax, threshold_rational;
#endif

class EigenSys
{
public :
  int order; Real_wp step;
  Real_wp dh_derivative, dh_second;
  Globatto<Real_wp> lob;
  
  Matrix<Complex_wp> rigid, Vh;
  Matrix<Real_wp, Symmetric, RowSymPacked> stiff; VectReal_wp mass;
  VectComplex_wp Yh, vec_tmp;
  
  EigenSys()
  {
    dh_derivative = 1e-9;
    dh_second = 1e-12;
  }
  
  void InitMatrix();
  void ComputeMatrix(const Real_wp& h);
  void EvaluateCflMin(const Real_wp& h, Real_wp& cfl, VectComplex_wp& eigen_vec);
  void EvaluateFunction(const VectReal_wp& h, VectReal_wp& der_cfl);
  void EvaluateJacobian(const VectReal_wp& h, Matrix<Real_wp>& jac_cfl);
  void FindInitGuess(VectReal_wp& h);
};

void EigenSys::InitMatrix()
{
  lob.ConstructQuadrature(order, lob.QUADRATURE_LOBATTO);
  lob.ComputeGradPhi();
  
  // matrice de masse
  mass.Reallocate(order);
  mass(0) = lob.Weights(0)*Real_wp(2);
  for (int i = 1; i < order; i++)
    mass(i) = lob.Weights(i);  
  
  // matrice de rigidite
  rigid.Reallocate(order, order);
  stiff.Reallocate(order+1, order+1); stiff.Fill(Real_wp(0));
  for (int i = 0; i <= order; i++)
    for (int j = i; j <= order; j++)
      for (int k = 0; k <= order; k++)
	stiff(i,j) += lob.Weights(k)*lob.GradPhi(i,k)*lob.GradPhi(j,k);
  
  Yh.Reallocate(order);
  Vh.Reallocate(order,order);
}

void EigenSys::ComputeMatrix(const Real_wp& h)
{
  // DISP(mass); DISP(stiff);
  rigid(0,0) = Complex_wp(stiff(0,0)*Real_wp(2) + Real_wp(2)*cos(h)*stiff(0,order), Real_wp(0));
  for (int i = 1; i < order; i++)
    {
      rigid(0,i) = stiff(0,i) + exp(-Iwp*h)*stiff(order,i);
      rigid(i,0) = conj(rigid(0,i));
      for (int j = 1; j < order; j++)
	rigid(i,j) = Complex_wp(stiff(i,j), Real_wp(0));
      
    }
  // DISP(rigid);

  Real_wp coef, c1, c2;
  // on met a l'echelle avec masse
  for (int i = 0; i < order; i++)
    for (int j = 0; j < order; j++)
      {
	c1 = sqrt(mass(i)); c2 = sqrt(mass(j));
	coef = c1*c2;
	rigid(i,j) /= coef;
      }
}

void EigenSys::EvaluateCflMin(const Real_wp& h, Real_wp& cfl, VectComplex_wp& eigen_vec)
{
  ComputeMatrix(h);
  GetEigenvaluesEigenvectors(rigid, Yh, Vh);
  // DISP(Yh);
  cfl = Real_wp(1);
  eigen_vec.Reallocate(order+1);
  for (int i = 0; i < order; i++)
    {
      Real_wp cfl_eval = Real_wp(Real_wp(2)/sqrt(abs(Yh(i))));
      if (cfl_eval < cfl)
	{
	  cfl = cfl_eval;
	  for (int j = 0; j < order; j++)
	    eigen_vec(j) = Vh(j,i);
	  
	  eigen_vec(order) = exp(Iwp*h)*eigen_vec(0);
	}
    }
}

void EigenSys::EvaluateFunction(const VectReal_wp& h, VectReal_wp& der_cfl)
{
  // derivee numerique
  Real_wp cfl_p, cfl_m;
  EvaluateCflMin(Real_wp(h(0)+dh_derivative), cfl_p, vec_tmp);
  EvaluateCflMin(Real_wp(h(0)-dh_derivative), cfl_m, vec_tmp);
  der_cfl(0) = (cfl_p-cfl_m)/(Real_wp(2)*dh_derivative);
}

void EigenSys::EvaluateJacobian(const VectReal_wp& h, Matrix<Real_wp>& der_second)
{
  // derivee numerique
  Real_wp cfl, cfl_p, cfl_m;
  EvaluateCflMin(Real_wp(h(0)+dh_second), cfl_p, vec_tmp);
  EvaluateCflMin(Real_wp(h(0)-dh_second), cfl_m, vec_tmp);
  EvaluateCflMin(Real_wp(h(0)), cfl, vec_tmp);
  der_second(0,0) = (cfl_p+cfl_m-Real_wp(2)*cfl)/(dh_second*dh_second);
}

void EigenSys::FindInitGuess(VectReal_wp& h)
{
  h(0) = step;
}

int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);
  
  VectComplex_wp Yh;
  // int order = atoi(argv[1]);
  int N = atoi(argv[1]);
  VectReal_wp cfl1D(N);
  ofstream file_out("cfl.dat"); file_out.precision(15);
  ofstream file_out_h("dx.dat"); file_out_h.precision(15);
  ofstream file_out_poly("poly.dat"); file_out_poly.precision(15);
  int nbx = 201;
  Real_wp deltax = Real_wp(1)/Real_wp(nbx-1);
  EigenSys sys;
  for (int order = 1; order <= N; order++)
    {
      sys.order = order;
      sys.InitMatrix();
      Real_wp cfl_min; VectComplex_wp eigen_vec;
      // if (order <= 20)
      // {
      // CFL atteinte pour h = 0 ou h = pi, suivant la parite de l'ordre
      if (order%2 == 0)
	sys.EvaluateCflMin(Real_wp(0), cfl_min, eigen_vec);
      else
	sys.EvaluateCflMin(pi_wp, cfl_min, eigen_vec);
      
      DISP(N); DISP(eigen_vec);
      for (int i = 0; i < nbx; i++)
	{
	  Real_wp xchap = deltax*i;
	  Complex_wp value(0);
	  for (int j = 0; j <= order; j++)
	    value += sys.lob.EvaluatePhi(j, xchap)*eigen_vec(j);
	  
	  file_out_poly<<real(value)<<"  ";
	}
      
      file_out_poly<<endl;
      
      // cfl_min = min(cfl_min, cfl);
	  
      // on fait une minimisation pour avoir la CFL exacte
      /* VectReal_wp h(1), dcfl(1), scale(1), rcontrol(20); 
	 IVect control(20);
	 MatrixFullReal_wp d2cfl(1,1);
	 cfl_min = Real_wp(1);
	 for (int i = 0; i <= (2*order); i++)
	 {
	      rcontrol.Fill(Real_wp(0)); control.Fill(0);
	      control(0) = 100;
	      control(3) = 1; control(4) = 0;
	      rcontrol(0) = 1e-6; rcontrol(1) = Real_wp(1);
	      // on cherche toutes les solutions dans l'intervalle [0,2pi]
	      sys.step = Real_wp(i)*pi_wp/order;
              sys.FindInitGuess(h);
	      int err = SolveMinpack(sys, h, dcfl, d2cfl, scale, control, rcontrol);
	      
	      sys.EvaluateCflMin(h(0), cfl);
	      cfl_min = min(cfl_min, cfl);
	      } */
      // }
      // else
      // {
      // on considere que le minimum n'est pas "loin" de la valeur obtenue pour h = 0
      // sys.EvaluateCflMin(Real_wp(0), cfl_min);
      // }
      DISP(cfl_min);
      
      // DISP(coef_cfl);
      cfl1D(order-1) = cfl_min;
      double cfld = toDouble(cfl_min);
      double dx = toDouble(sys.lob.Points(1)-sys.lob.Points(0));
      file_out<<cfld<<endl;
      file_out_h<<dx<<endl;
    }
  file_out.close();
  file_out_h.close();
  file_out_poly.close();
  
  return FinalizeMontjoie();
}
