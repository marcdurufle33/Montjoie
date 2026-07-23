#define MONTJOIE_WITH_TWO_DIM

#include "Quadrature/MontjoieQuadrature.hxx"
#include "Solver/MontjoieSolver.hxx"

using namespace Montjoie;

#ifdef MONTJOIE_WITH_MPFR
Real_wp kmax(0,MONTJOIE_DEFAULT_PRECISION), threshold_rational(0,MONTJOIE_DEFAULT_PRECISION);
#else
Real_wp kmax, threshold_rational;
#endif

class EigenSys : public NonLinearEquations_Base<Real_wp>
{
public :
  int order, type, nb_elem, elem_centre; Real_wp step;
  Real_wp dh_derivative, dh_second;
  Globatto<Real_wp> lob;
  
  Matrix<Complex_wp> rigid;
  Matrix<Real_wp, Symmetric, RowSymPacked> stiff; VectReal_wp mass;
  VectComplex_wp Yh; Matrix<Real_wp> ValExt;
  
  EigenSys()
  {
    dh_derivative = 1e-9;
    dh_second = 1e-12;
    nb_elem = 5; elem_centre = 2;
    //type = lob.QUADRATURE_LOBATTO;
    type = lob.QUADRATURE_GAUSS;
  }
  
  void InitMatrix();
  void ComputeMatrix(const Real_wp& h);
  void EvaluateCflMin(const Real_wp& h, Real_wp& cfl);
  void EvaluateFunction(const VectReal_wp& h, VectReal_wp& der_cfl);
  void EvaluateJacobian(const VectReal_wp& h, Matrix<Real_wp>& jac_cfl);
  void FindInitGuess(VectReal_wp& h);
};

void EigenSys::InitMatrix()
{
  lob.ConstructQuadrature(order, type);
  lob.ComputeGradPhi();
  
  ValExt.Reallocate(2, order+1);
  mass.Reallocate(order+1);
  for (int i = 0; i <= order; i++)
    {
      mass(i) = lob.Weights(i);
      ValExt(0,i) = lob.EvaluatePhi(i, Real_wp(0));
      ValExt(1,i) = lob.EvaluatePhi(i, Real_wp(1));
    }
  // DISP(ValExt);
  
  Matrix<Real_wp> stiffA;
  stiffA.Reallocate(nb_elem*(order+1), nb_elem*(order+1)); stiffA.Fill(Real_wp(0));
  // matrice de rigidite u avec v
  for (int num_elem = 0; num_elem < nb_elem; num_elem++)
    {
      // termes de rigidite
      for (int i = 0; i <= order; i++)
	for (int j = 0; j <= order; j++)
	  {
	    stiffA(i+num_elem*(order+1), j+num_elem*(order+1)) = 
	      lob.Weights(j)*lob.GradPhi(i,j) + ValExt(0,i)*ValExt(0,j)/2 - ValExt(1,i)*ValExt(1,j)/2;
	  }
      
      // termes de flux
      for (int i = 0; i <= order; i++)
	for (int j = 0; j <= order; j++)
	  {
	    // flux a gauche
	    if (num_elem > 0)
	      {
		stiffA(i+num_elem*(order+1), j+(num_elem-1)*(order+1)) = ValExt(0,i)*ValExt(1,j)/2;
	      }
	    else
	      {
		stiffA(i+num_elem*(order+1), j+(num_elem)*(order+1)) += ValExt(0,i)*ValExt(0,j)/2;
	      }
	    // flux a droite
	    if (num_elem < (nb_elem-1))
	      {
		stiffA(i+num_elem*(order+1), j+(num_elem+1)*(order+1)) = -ValExt(1,i)*ValExt(0,j)/2;
	      }
	    else
	      {
		stiffA(i+num_elem*(order+1), j+(num_elem)*(order+1)) += -ValExt(1,i)*ValExt(1,j)/2;
	      }
	  }
    }
  // DISP(stiffA);
  
  // computing R*R^t
  stiff.Reallocate(nb_elem*(order+1), nb_elem*(order+1));
  for (int i = 0; i < nb_elem*(order+1); i++)
    for (int j = i; j < nb_elem*(order+1); j++)
      {
	Real_wp vloc(0);
	for (int k = 0; k < nb_elem*(order+1); k++)
	  vloc = vloc + stiffA(i,k)*stiffA(j,k)/lob.Weights(k%(order+1));
	
	stiff(i,j) = vloc;
      }
  
  rigid.Reallocate(order+1, order+1);
}

void EigenSys::ComputeMatrix(const Real_wp& h)
{  
  // computing rigid, by periodicity
  rigid.Fill(Complex_wp(0,0));
  for (int i = 0; i <= order; i++)
    {
      for (int j = 0; j < nb_elem*(order+1); j++)
	{
	  int num_elem = j/(order+1);
	  int j1 = j%(order+1);
	  Complex_wp phase = exp(Iwp*Real_wp(num_elem-elem_centre)*h);
	  rigid(i,j1) += stiff(elem_centre*(order+1)+i,j)*phase;
	}
    }
  // DISP(rigid);
  
  Real_wp coef, c1, c2;
  // on met a l'echelle avec masse
  for (int i = 0; i <= order; i++)
    for (int j = 0; j <= order; j++)
      {
	c1 = sqrt(mass(i)); c2 = sqrt(mass(j));
	coef = c1*c2;
	rigid(i,j) /= coef;
      }
}

void EigenSys::EvaluateCflMin(const Real_wp& h, Real_wp& cfl)
{
  ComputeMatrix(h);
  GetEigenvalues(rigid, Yh);
  // DISP(Yh);
  cfl = Real_wp(1);
  for (int i = 0; i < order; i++)
    cfl = min(cfl, Real_wp(Real_wp(2)/sqrt(abs(Yh(i)))));
}

void EigenSys::EvaluateFunction(const VectReal_wp& h, VectReal_wp& der_cfl)
{
  // derivee numerique
  Real_wp cfl_p, cfl_m;
  EvaluateCflMin(Real_wp(h(0)+dh_derivative), cfl_p);
  EvaluateCflMin(Real_wp(h(0)-dh_derivative), cfl_m);
  der_cfl(0) = (cfl_p-cfl_m)/(Real_wp(2)*dh_derivative);
}

void EigenSys::EvaluateJacobian(const VectReal_wp& h, Matrix<Real_wp>& der_second)
{
  // derivee numerique
  Real_wp cfl, cfl_p, cfl_m;
  EvaluateCflMin(Real_wp(h(0)+dh_second), cfl_p);
  EvaluateCflMin(Real_wp(h(0)-dh_second), cfl_m);
  EvaluateCflMin(Real_wp(h(0)), cfl);
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
  VectReal_wp cfl1D(N+1); cfl1D(0) = Real_wp(1);
  ofstream file_out("cfl.dat"); file_out.precision(7);
  ofstream file_out_h("dx.dat"); file_out_h.precision(7);
  EigenSys sys;
  for (int order = 1; order <= N; order++)
    {
      sys.order = order;
      sys.InitMatrix();
      Real_wp cfl, cfl_min;
      // if (order <= 20)
      // {
      // CFL atteinte pour h = 0 ou h = pi, suivant la parite de l'ordre
      // if (order%2 == 0)
      Real_wp h = 1e-5;
      sys.ComputeMatrix(h);
      GetEigenvalues(sys.rigid, Yh);
      Real_wp err(1);
      for (int i = 0; i < Yh.GetM(); i++)
	err = min(err, Real_wp(abs(Real_wp(1)-abs(Yh(i))/(h*h))));
      
      /// DISP(Yh); DISP(err);
      if (order%2 == 0)
	{
	  Real_wp alpha = pow(h, Real_wp(2*(order+1)))/err;
	  // DISP(alpha);
	  alpha = pow(h, Real_wp(2*order))/err;
	  // DISP(alpha);
	}
      else
	{
	  Real_wp alpha = pow(h, Real_wp(2*order))/err;
	  // DISP(alpha);
	}
      
      sys.EvaluateCflMin(Real_wp(0), cfl_min);
      // else
      // sys.EvaluateCflMin(1.01*pi_wp, cfl); // DISP(cfl);
      // sys.EvaluateCflMin(2*pi_wp, cfl); // DISP(cfl);
      sys.EvaluateCflMin(pi_wp, cfl); // DISP(cfl);
      
      cfl_min = min(cfl_min, cfl);
	  
      // on fait une minimisation pour avoir la CFL exacte
      VectReal_wp h_(1), dcfl(1), scale(1), rcontrol(20); 
      IVect control(20);
      Matrix<Real_wp> d2cfl(1,1);
      cfl_min = Real_wp(1);
      for (int i = 0; i <= 1000; i++)
	{
	  rcontrol.Fill(Real_wp(0)); control.Fill(0);
	  control(0) = 1000;
	  control(3) = 1; control(4) = 0;
	  rcontrol(0) = 1e-10; rcontrol(1) = Real_wp(1);
	  // on cherche toutes les solutions dans l'intervalle [0,2pi]
	  sys.step = Real_wp(i)*pi_wp/1000;
          sys.FindInitGuess(h_);
	  int err = SolveMinpack(sys, h_, dcfl, d2cfl, scale, control, rcontrol);
	  
	  sys.EvaluateCflMin(h_(0), cfl); // DISP(h_(0)); DISP(cfl);
	  cfl_min = min(cfl_min, cfl);
	}
      DISP(cfl_min);
      
      // DISP(coef_cfl);
      cfl1D(order) = cfl_min;
      double cfld = toDouble(cfl_min);
      double dx(1);
      if (order > 0) 
	dx = toDouble(sys.lob.Points(1)-sys.lob.Points(0));
      file_out<<cfld<<endl;
      file_out_h<<dx<<endl;
    }
  file_out.close();
  file_out_h.close();
  
  return FinalizeMontjoie();
}
