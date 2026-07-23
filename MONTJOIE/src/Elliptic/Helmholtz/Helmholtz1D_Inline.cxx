#ifndef MONTJOIE_FILE_HELMHOLTZ_1D_INLINE_CXX

namespace Montjoie
{

  inline bool HelmholtzEquation1D::ComputeDFjm1()
  {
    return false; 
  }
  

  inline bool LaplaceEquation1D::ComputeDFjm1()
  {
    return false; 
  }


  /*******************
   * VarHelmholtz_1D *
   *******************/
  
  
  template<class Complexe>
  inline void VarHelmholtz_1D<Complexe>::GetIkwave(Real_wp& ikwave)
  {
    ikwave = 1.0;
  }
  
  
  template<class Complexe>
  inline void VarHelmholtz_1D<Complexe>::GetIkwave(Complex_wp& ikwave)
  {
    ikwave = Iwp*kwave;
  }
  

  template<class Complexe>    
  inline void VarHelmholtz_1D<Complexe>
  ::ComputeElementaryMatrix(int i, IVect& num, VirtualMatrix<Real_wp>& A,
                            CondensationBlockSolver_Base<Real_wp>&,
			    const GlobalGenericMatrix<Real_wp>& nat_mat)
  {
    abort();
  }
  
  
  template<class Complexe>    
  inline void VarHelmholtz_1D<Complexe>
  ::ComputeElementaryMatrix(int i, IVect& num, VirtualMatrix<Complex_wp>& A,
                            CondensationBlockSolver_Base<Complex_wp>&,
			    const GlobalGenericMatrix<Complex_wp>& nat_mat)
  {
    ComputeElementaryMatrixGen(i, num, A, nat_mat);
  }


  inline GaussianSource<Dimension1>::GaussianSource()
  {
    x0 = 0.0;
    radius = 1.0;
    alpha = 1.0;
    beta = 1.0;
  }

  
  inline Real_wp GaussianSource<Dimension1>::GetAmplitude(const Real_wp& x) const
  {
    Real_wp aj = beta*exp(-alpha*square(x-x0));
    return aj;
  }
  
  
  template<class T>
  inline void GaussianSourceField<T, Dimension1>
  ::EvaluateFunction(const Real_wp& x, T& f) const
  {
    Real_wp aj = this->GetAmplitude(x);
    f = aj;
  }
  
  
  template<class T>
  inline void GaussianSourceField<T, Dimension1>
  ::EvaluateFunctionGradient(const Real_wp& x, T& f, TinyVector<T, 1>& grad_f) const
  {
    Real_wp aj = this->GetAmplitude(x);
    f = aj;
    grad_f(0) = -2.0*this->alpha*(x-this->x0)*aj; 
  }

  
  inline void EllipticProblem<LaplaceEquation1D>
  ::ComputeElementaryMatrix(int i, IVect& num, VirtualMatrix<Real_wp>& A,
                            CondensationBlockSolver_Base<Real_wp>&,
			    const GlobalGenericMatrix<Real_wp>& nat_mat)
  {
    this->ComputeElementaryMatrixGen(i, num, A, nat_mat);
  }


  inline void EllipticProblem<LaplaceEquation1D>
  ::ComputeElementaryMatrix(int i, IVect& num, VirtualMatrix<Complex_wp>& A,
                            CondensationBlockSolver_Base<Complex_wp>&,
			    const GlobalGenericMatrix<Complex_wp>& nat_mat)
  {
    this->ComputeElementaryMatrixGen(i, num, A, nat_mat);
  }
  
}

#define MONTJOIE_FILE_HELMHOLTZ_1D_INLINE_CXX
#endif
