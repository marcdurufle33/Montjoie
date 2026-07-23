#ifndef MONTJOIE_FILE_IMPEDANCE_HELMHOLTZ_INLINE_CXX

namespace Montjoie
{
  
  template<class Dimension>
  inline ImpedanceABC<Complex_wp, HelmholtzEquation<Dimension> >
  ::ImpedanceABC(const EllipticProblem<HelmholtzEquation<Dimension> >& var)
    : ImpedanceABC_Helm<Dimension>(var) 
  {
  }

  
  template<class Dimension>
  inline ImpedanceABC<Complex_wp, HelmholtzEquationDG<Dimension> >
  ::ImpedanceABC(const EllipticProblem<HelmholtzEquationDG<Dimension> >& var)
    : ImpedanceABC_Helm<Dimension>(var) 
  {
  }


  //! Constructor with the given Helmholtz problem
  template<class Dimension> template<class TypeEquation>
  inline ImpedanceABC_Helm<Dimension>
  ::ImpedanceABC_Helm(const EllipticProblem<TypeEquation>& var)
    : ImpedanceFunction_Base<Complex_wp, Dimension>(var),
      var_helm(var), var_boundary(var), var_problem(var)
    , var_gibc(var.var_gibc)
  {
    InitDefaultValues();
  }
  

  template<class Dimension>
  inline bool ImpedanceABC_Helm<Dimension>::PresenceGradient() const
  {
    if ((var_boundary.GetOrderAbsorbingCondition() == 1)
	&& (!var_boundary.grazing_abc))
      return false;
    
    return true;
  }
  
#ifdef MONTJOIE_WITH_HIGH_CONDUCTIVITY_MODEL
  template<class Dimension>
  inline bool ImpedanceHighConductivity<Complex_wp, HelmholtzEquation<Dimension> >::PresenceGradient() const
  {
    return true;
  }
#endif
  
  //! default constructor
  template<class Dimension> template<class TypeEquation>
  inline VarGeneralizedImpedance_Helm<Dimension>
  ::VarGeneralizedImpedance_Helm(EllipticProblem<TypeEquation>& var)
    : VarGeneralizedImpedance_Base(var), var_problem(var), var_helm(var), var_boundary(var)
  {
  }
  

  template<class Dimension>
  inline void VarGeneralizedImpedance_Helm<Dimension>
  ::AddGibcTerms(const Real_wp& alpha, const GlobalGenericMatrix<Real_wp>& nat_mat,
                 VirtualMatrix<Real_wp>& mat_sp, int offset_row, int offset_col)
  {
    abort(); 
  }

  
  template<class Dimension>
  inline VarGeneralizedImpedance<HelmholtzEquation<Dimension> > 
  ::VarGeneralizedImpedance(EllipticProblem<HelmholtzEquation<Dimension> >& var)
    : VarGeneralizedImpedance_Helm<Dimension>(var) 
  {
  }


  template<class Dimension>
  inline VarGeneralizedImpedance<HelmholtzEquationDG<Dimension> > 
  ::VarGeneralizedImpedance(EllipticProblem<HelmholtzEquationDG<Dimension> >& var)
    : VarGeneralizedImpedance_Helm<Dimension>(var) 
  {
  }

}

#define MONTJOIE_FILE_IMPEDANCE_HELMHOLTZ_INLINE_CXX
#endif
