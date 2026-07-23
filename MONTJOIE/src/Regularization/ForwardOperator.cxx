#ifndef ITREG_FILE_FORWARD_OPERATOR_CXX

namespace itreg
{

  //! Evaluates jacobian matrix
  template<class T, class VectorRhs, class VectorSol>
  void ForwardOperator<T, VectorRhs, VectorSol>
  ::FullDerivative(const VectorSol& x, const VectorRhs& y, Matrix<T>& DF)
  {
    VectorSol Ones(x);
    VectorRhs F_ones(y);
    Ones.Zero();

    DF.Reallocate(y.GetM(), x.GetM());
    // computing each column of DF by using Derivative
    for (int j = 0; j < x.GetM(); j++)
      {
	Ones(j) = T(1);
	Derivative(Ones, F_ones);
	
	for (int i = 0; i < y.GetM(); i++)
	  DF(i, j) = F_ones(i);
	
	Ones(j) = T(0);
      }
  }


  //! Evaluates jacobian matrix DF(x) and F(x)
  template<class T, class VectorRhs, class VectorSol>
  void ForwardOperator<T, VectorRhs, VectorSol>
  ::EvaluateJacobian(const VectorSol& x, VectorRhs& y, Matrix<T>& DF)
  {
    // using FullDerivative
    this->Init(x);
    this->Evaluate(y);
    this->FullDerivative(x, y, DF);    
  }

  
  //! Checks if the derivative and ajoint derivative are correct
  template<class T, class VectorRhs, class VectorSol>
  void ForwardOperator<T, VectorRhs, VectorSol>::CheckDerivatives(const Vector<T>& c0, const T& h)
  {
    Vector<T> f, f0; Matrix<T> DF, DF_full;
    // testing Derivative by comparing with second-order finite difference scheme for DF
    cout << "Testing DF" << endl;
    for (int j = 0; j < c0.GetM(); j++)
      {
	DISP(j);
	Vector<T> ctmp(c0), fm, fp;
	ctmp(j) += h;
	
	this->Init(ctmp);
	this->Evaluate(fp);
	
	ctmp(j) -= 2.0*h;
	
	this->Init(ctmp);
	this->Evaluate(fm);
	
	Vector<T> DF_num;
	DF_num = (fp - fm) / (2.0*h);
	
	Vector<T> DF_exact, Ones;
	Ones.Reallocate(c0.GetM());
	Ones.Zero(); Ones(j) = 1.0;
	this->Init(c0);
	this->Derivative(Ones, DF_exact);
	if (j == 0)
	  DF.Reallocate(DF_exact.GetM(), c0.GetM());
	
	SetCol(DF_exact, j, DF);
	for (int i = 0; i < DF_exact.GetM(); i++)
	  if ((abs(DF_exact(i) - DF_num(i)) > h*sqrt(h)) || isnan(abs(DF_exact(i) - DF_num(i))))
	    {
	      DISP(j); DISP(i); DISP(DF_exact(i)); DISP(DF_num(i));
              DISP(DF_exact(i) - DF_num(i));
	      cout << "Derivative inexact" << endl;
	      abort();
	    }
      }

    // testing EvaluateJacobian
    cout << "Testing EvaluateJacobian" << endl;
    this->EvaluateJacobian(c0, f, DF_full);

    this->Init(c0);
    this->Evaluate(f0);

    for (int i = 0; i < f0.GetM(); i++)
      {
	if ((abs(f(i)-f0(i)) > h*sqrt(h)) || isnan(abs(f(i)-f0(i))))
	  {
	    DISP(i); DISP(f(i)); DISP(f0(i));
	    cout << "EvaluateJacobian inexact" << endl;
	    abort();
	  }

	for (int j = 0; j < DF.GetN(); j++)
	  if ((abs(DF(i, j) - DF_full(i, j)) > h*sqrt(h)) || isnan(abs(DF(i, j) - DF_full(i, j))))
	    {
	      DISP(i); DISP(j); DISP(DF(i, j)); DISP(DF_full(i, j));
	      cout << "EvaluateJacobian inexact" << endl;
	      abort();
	    }
      }

    // testing Adjoint
    cout << "Testing DF*" << endl;
    this->Init(c0);
    for (int j = 0; j < DF.GetM(); j++)
      {
	DISP(j);
	Vector<T> DFT_exact, Ones;
	Ones.Reallocate(DF.GetM());
	Ones.Zero(); Ones(j) = 1.0;
	this->Adjoint(Ones, DFT_exact);
	
	for (int i = 0; i < c0.GetM(); i++)
	  if ((abs(DFT_exact(i) - DF(j, i)) > h*sqrt(h)) || isnan(abs(DFT_exact(i) - DF(j, i))))
	    {
	      DISP(i); DISP(j); DISP(DFT_exact(i)); DISP(DF(j, i));
	      DISP(DFT_exact);
	      cout << "Adjoint inexact" << endl;
	      abort();
	    }
      }
  }
  
}

#define ITREG_FILE_FORWARD_OPERATOR_CXX
#endif

