#include "Share/MontjoieCommon.hxx"
#include "Solver/MontjoieSolver.hxx"

using namespace Montjoie;

class FittingExample : public VirtualLeastSquaresFunction<Real_wp>
{
private:
  VectReal_wp pts_ti, yi;
  
public:
  FittingExample()
  {
    this->m_ = 10;
    this->n_ = 4;
    
    Linspace(Real_wp(0), Real_wp(1), this->m_, pts_ti);
    
    Real_wp x1 = 4, x2 = 5, x3 = 4, x4 = 4;
    yi.Reallocate(this->m_);
    for (int i = 0; i < this->m_; i++)
      {
	Real_wp t = pts_ti(i);
	yi(i) = x3*exp(x1*t) + x4*exp(x2*t);
      }
    
    DISP(pts_ti); DISP(yi);
  }

  void FindInitGuess(Vector<Real_wp>& x)
  {
    x.Reallocate(this->n_);
    x(0) = -1;
    x(1) = -2;
    x(2) = 1;
    x(3) = -1;
  }

  void EvaluateF(const Vector<Real_wp>& x, Vector<Real_wp>& f)
  {
    for (int i = 0; i < this->m_; i++)
      {
	Real_wp t = pts_ti(i);
	f(i) = yi(i) - x(2)*exp(x(0)*t) - x(3)*exp(x(1)*t);
      }
  }

  void EvaluateJacobian(const Vector<Real_wp>& x, Vector<Real_wp>& f,
			Matrix<Real_wp, General, ColMajor>& fjac)
  {
    for (int i = 0; i < this->m_; i++)
      {
	Real_wp t = pts_ti(i);
	f(i) = yi(i) - x(2)*exp(x(0)*t) - x(3)*exp(x(1)*t);
	fjac(i, 0) = -t*x(2)*exp(x(0)*t);
	fjac(i, 1) = -t*x(3)*exp(x(1)*t);
	fjac(i, 2) = -exp(x(0)*t);
	fjac(i, 3) = -exp(x(1)*t);
      }    
  }
  
};


class FittingExampleTiny
  : public VirtualLeastSquaresFunction<Real_wp, TinyVector<Real_wp, 4>, TinyVector<Real_wp, 10>,
				       TinyMatrix<Real_wp, General, 10, 4>,
				       TinyMatrix<Real_wp, Symmetric, 4, 4> >
{
private:
  VectReal_wp pts_ti, yi;
  
public:
  FittingExampleTiny()
  {
    this->m_ = 10;
    this->n_ = 4;
    
    Linspace(Real_wp(0), Real_wp(1), this->m_, pts_ti);
    
    Real_wp x1 = 4, x2 = 5, x3 = 4, x4 = 4;
    yi.Reallocate(this->m_);
    for (int i = 0; i < this->m_; i++)
      {
	Real_wp t = pts_ti(i);
	yi(i) = x3*exp(x1*t) + x4*exp(x2*t);
      }
    
    DISP(pts_ti); DISP(yi);
  }

  void FindInitGuess(TinyVector<Real_wp, 4>& x)
  {
    x(0) = -1;
    x(1) = -2;
    x(2) = 1;
    x(3) = -1;
  }

  void EvaluateF(const TinyVector<Real_wp, 4>& x, TinyVector<Real_wp, 10>& f)
  {
    for (int i = 0; i < this->m_; i++)
      {
	Real_wp t = pts_ti(i);
	f(i) = yi(i) - x(2)*exp(x(0)*t) - x(3)*exp(x(1)*t);
      }
  }

  void EvaluateJacobian(const TinyVector<Real_wp, 4>& x, TinyVector<Real_wp, 10>& f,
			TinyMatrix<Real_wp, General, 10, 4>& fjac)
  {
    for (int i = 0; i < this->m_; i++)
      {
	Real_wp t = pts_ti(i);
	f(i) = yi(i) - x(2)*exp(x(0)*t) - x(3)*exp(x(1)*t);
	fjac(i, 0) = -t*x(2)*exp(x(0)*t);
	fjac(i, 1) = -t*x(3)*exp(x(1)*t);
	fjac(i, 2) = -exp(x(0)*t);
	fjac(i, 3) = -exp(x(1)*t);
      }    
  }
  
};

int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);
  
  FittingExample fct;
  Vector<Real_wp> xsol;
  SolveLeastSquaresLvm(fct, xsol, 1e-13, 1e-13, 1000, 2, 1e-3);
  
  Real_wp threshold = 1e-11;
  DISP(xsol);
  if ((abs(xsol(0) - Real_wp(5)) > threshold) || (abs(xsol(1) - Real_wp(4)) > threshold)
      || (abs(xsol(2) - Real_wp(4)) > threshold) || (abs(xsol(3) - Real_wp(4)) > threshold))
    {
      cout << "SolveLeastSquaresLvm incorrect" << endl;
      abort();
    }


  FittingExampleTiny fct_tiny;
  TinyVector<Real_wp, 4> xsol_tiny;
  SolveLeastSquaresLvm(fct_tiny, xsol_tiny, 1e-13, 1e-13, 1000, 2, 1e-3);
  
  DISP(xsol_tiny);
  if ((abs(xsol_tiny(0) - Real_wp(5)) > threshold) || (abs(xsol_tiny(1) - Real_wp(4)) > threshold)
      || (abs(xsol_tiny(2) - Real_wp(4)) > threshold) || (abs(xsol_tiny(3) - Real_wp(4)) > threshold))
    {
      cout << "SolveLeastSquaresLvm incorrect" << endl;
      abort();
    }

  cout << "All tests passed successfully" << endl;
  
  return FinalizeMontjoie();
}
