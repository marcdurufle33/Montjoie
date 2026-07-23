#ifndef MONTJOIE_FILE_ADVECTION_PROBLEM_HXX

#include "Algebra/TinyBandMatrix.cxx"

namespace Montjoie
{
  //! base class for advection equation 
  /*!
    We are considering advection equation du/dt + c du/dx + 1/2 c' u = 0
    where c' is the derivative of the velocity c 
   */
  class AdvectionEquation : public GenericWaterWaveEquation
  {
  public :
    // evaluation of c and c'
    VectReal_wp Cquad, CxQuad;
    
    AdvectionEquation();
    
    template<class GenericPb>
    void ComputePhysicalProperty(const GenericPb& var, double xmin, double xmax);

  };
  
  
  class AdvectionEquationDG;
  
  //! operator for the LDG stiffness matrix for advection equation
  class AdvectionStiffnessOperatorDG
  {
  protected :
    AdvectionEquationDG& var;
    
  public :    
    AdvectionStiffnessOperatorDG(AdvectionEquationDG& var_kdv);
    
    template<class Vector1, class Vector2>
    void EvaluateNumericalFlux(int i, int k, const Vector1& Um,
			       const Vector1& Up, Vector2& flux) const;
    
    template<class Vector1, class Vector2, class Matrix1>
    void EvaluateDerivativeNumericalFlux(int i, int k, const Vector1& Um,
					 const Vector1& Up, Vector2& flux,
					 Matrix1& df_dum, Matrix1& df_dup) const;
    
    template<class Vector1, class Vector2>
    void EvaluateFunction(int i, const Vector1& u, const Vector1& du,
			  Vector2& g, Vector2& h) const;
    
    template<class Vector1, class Vector2, class Matrix1>
    void EvaluateDerivative(int i, const Vector1& u, const Vector1& du,
                            Vector2& g, Vector2& h, Matrix1& df_du,
                            Matrix1& dg_du, Matrix1& df_ddu, Matrix1& dg_ddu) const;
    
  };
  
  
  //! class for LDG formulation of the advection equation
  class AdvectionEquationDG : public AdvectionEquation
  {
  public :
    //! penalization parameter
    Real_wp alpha_lax;
    IdentityOperatorDg1D<AdvectionEquationDG> var_mass;
    AdvectionStiffnessOperatorDG var_stiff;
    enum { nb_unknowns = 1, nb_time_unknowns = 1};
    
    AdvectionEquationDG();
    AdvectionEquationDG(const AdvectionEquationDG& var);
    
    static const bool FormulationDG = true;
    static bool LinearMassTerm() { return true; }
    
  };


  class AdvectionEquationFem;
  
  //! operator for the FEM stiffness matrix for advection equation  
  class AdvectionStiffnessOperatorFem : public IdentityOperatorDg1D<AdvectionEquationFem>
  {
  protected :
    AdvectionEquationFem& var;
    
  public :    
    AdvectionStiffnessOperatorFem(AdvectionEquationFem& var_kdv);
    
    template<class Vector1, class Vector2>
    void EvaluateFunction(int i, const Vector1& u, const Vector1& du,
			  Vector2& g, Vector2& h) const;
    
    template<class Vector1, class Vector2, class Matrix1>
    void EvaluateDerivative(int i, const Vector1& u, const Vector1& du,
                            Vector2& g, Vector2& h, Matrix1& df_du,
                            Matrix1& dg_du, Matrix1& df_ddu, Matrix1& dg_ddu) const;
    
  };
  
  
  //! class for the continuous formulation of the advection equation  
  class AdvectionEquationFem : public AdvectionEquation
  {
  public :
    IdentityOperatorDg1D<AdvectionEquationFem> var_mass;
    AdvectionStiffnessOperatorFem var_stiff;
    enum { nb_unknowns = 1, nb_time_unknowns = 1};
    
    static const bool FormulationDG = false;
    
    AdvectionEquationFem();
    AdvectionEquationFem(const AdvectionEquationFem& var);
    
    static bool LinearMassTerm() { return true; }
    
  };

  
  //! class for the resolution of advection equation with finite difference method
  template<>
  class FiniteDifferenceScheme1D<AdvectionEquation>
    : public FiniteDifferenceScheme1D_Base<AdvectionEquation>
  {
  protected :
    
  public :    
    
    FiniteDifferenceScheme1D();

    void AdvancePredictorScheme(Real_wp t, int iter_num)
    { abort(); }
    
    void AddFunction(const Real_wp& alpha, const VectReal_wp& Xn,
                     VectReal_wp& ProdXn, bool invert_mass = true)
    {
      abort(); 
    }

    void GetFunction(const VectReal_wp& Xn, VectReal_wp& ProdXn)
    { abort(); }
    
    void ComputeAndFactoriseJacobian(const Real_wp& alpha,
                                     const VectReal_wp& Un_demi)
    { abort(); }

    void SolveJacobian(VectReal_wp& rhs)
    { abort(); }

    void InitJacobian(int incx) { abort(); }

    void AddStiffJacobian(const Real_wp& alpha,
                          const VectReal_wp& Un_demi,
                          int m = 0, int n = 0, int incx = 1, int incy = 1)
    { abort(); }

    void FactoriseJacobian(int incx)
    { abort(); }
    
    void AddMass(int m = 0, int n = 0, int incx = 1, int incy = 1)
    { abort(); }

    void SolveGaussJacobian(VectReal_wp& rhs, int incx)
    { abort(); }

  };
  
  
  //! class for the resolution of the advection equation with spectral method
  class AdvectionEquationSpectral : public AdvectionEquation
  {
  public :
    enum { nb_unknowns = 1, nb_time_unknowns = 1};
    typedef AdvectionEquation TypeEquationDiff;
    
  protected :
    bool constant_advec;
    VectReal_wp Xreal, dXreal, ProdReal;
    VectComplex_wp XnHat, dXn;
    
  public :
    AdvectionEquationSpectral();
    
    void ComputePhysicalProperty(FftRealInterface& fft, const VectReal_wp& nu,
				 double xmin, double xmax, int N);
    
    void InitSinus(const Real_wp& dt, const VectReal_wp& nu) { abort(); }

    void ApplyMass(const VectReal_wp&,
                   FftRealInterface& fft, VectReal_wp& nu, VectReal_wp&) { abort(); }
    
    void AddFunction(const Real_wp&, const VectReal_wp&,
                     FftRealInterface& fft, VectReal_wp& nu,
                     VectReal_wp&, bool solve_mass = true) { abort(); }
    
    void EvaluateFunction(const Real_wp& t, const VectComplex_wp& Xn,
			  FftRealInterface& fft, VectReal_wp& nu,
			  VectComplex_wp& ProdXn);

    void EvaluateFunction(const Real_wp& t, const VectReal_wp& Xn,
			  FftRealInterface& fft, VectReal_wp& nu,
			  VectReal_wp& ProdXn);
    
  };

}

#define MONTJOIE_FILE_ADVECTION_PROBLEM_HXX
#endif
