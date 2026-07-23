#ifndef MONTJOIE_FILE_HARMONIC_MAXWELL_2D_INLINE_CXX

namespace Montjoie
{
   
  //! penalization coefficient
  template<class TypeEquation, class T0>
  inline void HarmonicMaxwellEquation_2D::
  GetPenalizationCoef(const EllipticProblem<TypeEquation>& var,
		      int i, const GlobalGenericMatrix<T0>& nat_mat, int ref, T0& coef)
  {
    coef = var.delta_penalization*nat_mat.GetCoefStiffness()
      /var.coefficient_impedance_absorbing(ref);
  }
  
  
  /**********************
   * HarmonicMaxwell_2D *
   **********************/
  
  
  //! constructor with a given problem
  template<class Complexe> template<class TypeEquation>
  inline HarmonicMaxwell_2D<Complexe>
  ::HarmonicMaxwell_2D(EllipticProblem<TypeEquation>& var)
    : var_problem(var), var_boundary(var)
  {
    SetComplexOne(epsilon0);
    SetComplexOne(mu0);
    SetComplexOne(invMu0);
  }
  
  
  //! returns mass coefficients involved in mass matrix
  template<class Complexe> template<class T0, class Prop>
  inline void HarmonicMaxwell_2D<Complexe>::
  GetMassMatrixCoef(int num_elem, int j,
		    const GlobalGenericMatrix<T0>& nat_mat, int m, int n,
		    TinyMatrix<T0, Prop, 2, 2>& Bmass) const
  {
    Bmass = this->Glob_matMass_Bh(num_elem)(j);
    Mlt(nat_mat.GetCoefMass(), Bmass);
  }
  
  
  //! returns mass coefficients involved in stiffness matrix
  template<class Complexe> template<class T0>
  inline void HarmonicMaxwell_2D<Complexe>::
  GetStiffMatrixCoef(int num_elem, int j,
		     const GlobalGenericMatrix<T0>& nat_mat, int m, int n, T0& Astiff) const
  {
    if (j < this->Glob_matMass_Dh(num_elem).GetM())
      {
	Astiff = this->Glob_matMass_Dh(num_elem)(j);
	if (var_problem.FirstOrderFormulation())
	  Astiff *= nat_mat.GetCoefMass();
	else
	  Astiff *= nat_mat.GetCoefStiffness();
      }
    else
      Astiff = 0.0;
  }
  
      
  /****************
   * ImpedanceABC *
   ****************/
  

  template<class TypeEquation>
  inline ImpedanceABC<Complex_wp, HarmonicMaxwellEquation_2D>
  ::ImpedanceABC(const EllipticProblem<TypeEquation>& var)
    : ImpedanceFunction_Base<Complex_wp, Dimension2>(var), var_problem(var), var_maxwell(var)
  {
  }


  inline FemMatrixFreeClass<Complex_wp, HarmonicMaxwellEquation_2D>
  ::FemMatrixFreeClass(const EllipticProblem<HarmonicMaxwellEquation_2D>& var_)
    : FemMatrixFreeClass_Eq<Complex_wp, HarmonicMaxwellEquation_2D>(var_)
  {
  }
          
} // end namespace

#define MONTJOIE_FILE_HARMONIC_MAXWELL_2D_INLINE_CXX
#endif
