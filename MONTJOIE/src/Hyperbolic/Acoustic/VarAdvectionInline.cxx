#ifndef MONTJOIE_FILE_VAR_ADVECTION_INLINE_CXX

namespace Montjoie
{

  template<class T, class Dim>
  inline bool AdvectionEquation_Base<T, Dim>::SymmetricGlobalMatrix()
  {
    return false; 
  }
  
  
  template<class T, class Dim>
  inline bool AdvectionEquation_Base<T, Dim>::SymmetricElementaryMatrix()
  {
    return false;
  }    
  
  
  template<class T, class Dim>
  inline bool AdvectionEquation_Base<T, Dim>::ComputeDFjm1()
  {
    return true;
  }

  
  template<class T, class Dim>
  inline const Vector<TinyMatrix<T, Symmetric,
				 Dim::dim_N, Dim::dim_N> >&
  AdvectionEquation_Base<T, Dim>
  ::GetStiffCoefficient(const VarHelmholtz_Base<T, Dimension>& var, int num_elem, int ref)
  {
    return var.Glob_matMass_Bh(num_elem); 
  }
  
  
  template<class T, class Dim>
  inline const Vector<T>& AdvectionEquation_Base<T, Dim>
  ::GetMassCoefficient(const VarHelmholtz_Base<T, Dimension>& var, int num_elem, int ref)
  { 
    return var.Glob_matMass_Dh(num_elem);
  }
  
  
  //! constructor with a given problem
  template<class Complexe, class Dimension> template<class TypeEquation>
  inline ImpedanceABC<Complexe, AdvectionEquation<Dimension> >
  ::ImpedanceABC(const EllipticProblem<TypeEquation>& var)
    : ImpedanceFunction_Base<Complexe, Dimension>(var), var_helm(var), var_problem(var) 
  {
  }
  
  
  //! constructor with a given problem
  template<class Dimension> template<class TypeEquation>
  inline ImpedanceABC<Complex_wp, HarmonicAdvectionEquation<Dimension> > 
  ::ImpedanceABC(const EllipticProblem<TypeEquation>& var)
    : ImpedanceFunction_Base<Complex_wp, Dimension>(var), var_helm(var), var_problem(var) 
  {
  }
  

  template<class Dimension>
  inline EllipticProblem<AdvectionEquation<Dimension> >::EllipticProblem()
  {
    this->dg_formulation = ElementReference_Base::CONTINUOUS; 
  }


  template<class Dimension>
  inline EllipticProblem<HarmonicAdvectionEquation<Dimension> >::EllipticProblem()
  {
    this->dg_formulation = ElementReference_Base::CONTINUOUS; 
  }
  

  template<class T, class Dim>
  inline bool DoubleAdvectionEquation_Base<T, Dim>::ComputeDFjm1()
  {
    return true; 
  }

  
  template<class T, class Dim>
  inline bool DoubleAdvectionEquation_Base<T, Dim>::SymmetricGlobalMatrix()
  {
    return false;
  }
    
  
  template<class T, class Dim>
  inline bool DoubleAdvectionEquation_Base<T, Dim>::SymmetricElementaryMatrix()
  {
    return false;
  }    
  
  
  template<class T, class Dim>
  inline const Vector<TinyMatrix<T, Symmetric, Dim::dim_N, Dim::dim_N> >&
  DoubleAdvectionEquation_Base<T, Dim>
  ::GetStiffCoefficient(const VarHelmholtz_Base<T, Dimension>& var, int num_elem, int ref)
  {
    return var.Glob_matMass_Bh(num_elem); 
  }
  
    
  template<class T, class Dim>
  inline const Vector<T>& DoubleAdvectionEquation_Base<T, Dim>::
  GetMassCoefficient(const VarHelmholtz_Base<T, Dimension>& var, int num_elem, int ref)
  { 
    return var.Glob_matMass_Dh(num_elem); 
  }
  
  
  template<class Complexe, class Dimension> template<class TypeEquation>
  inline ImpedanceABC_DoubleAdvec<Complexe, Dimension>
  ::ImpedanceABC_DoubleAdvec(const EllipticProblem<TypeEquation>& var)
    : ImpedanceFunction_Base<Complexe, Dimension>(var), var_helm(var), var_problem(var) 
  {
  }
  

  template<class Complexe, class Dimension> template<class TypeEquation>
  inline ImpedanceABC<Complexe, DoubleAdvectionEquation<Dimension> >
  ::ImpedanceABC(const EllipticProblem<TypeEquation>& var)
    : ImpedanceABC_DoubleAdvec<Complexe, Dimension>(var) 
  {
  }


  template<class Dimension> template<class TypeEquation>
  inline ImpedanceABC<Complex_wp, HarmonicDoubleAdvectionEquation<Dimension> > 
  ::ImpedanceABC(const EllipticProblem<TypeEquation>& var)
    : ImpedanceABC_DoubleAdvec<Complex_wp, Dimension>(var)
  {
  }

  
  template<class Dimension>
  inline EllipticProblem<DoubleAdvectionEquation<Dimension> >::EllipticProblem()
  {
    this->dg_formulation = ElementReference_Base::CONTINUOUS; 
  }

  
  template<class Dimension>
  inline EllipticProblem<HarmonicDoubleAdvectionEquation<Dimension> >::EllipticProblem()
  {
    this->dg_formulation = ElementReference_Base::CONTINUOUS;
  }

}

#define MONTJOIE_FILE_VAR_ADVECTION_INLINE_CXX
#endif
