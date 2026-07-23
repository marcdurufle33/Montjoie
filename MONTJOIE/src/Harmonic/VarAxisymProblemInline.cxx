#ifndef MONTJOIE_FILE_VAR_AXISYM_PROBLEM_INLINE_CXX

namespace Montjoie
{
  
#ifdef MONTJOIE_WITH_TWO_DIM
  /********************
   * VarAxisymProblem *
   ********************/

  
  //! constructor with a given problem
  template<class TypeEquation>  
  inline VarAxisymProblem::VarAxisymProblem(EllipticProblem<TypeEquation>& var)
    : var_problem(var), var_boundary(var), var_source(var)
    , var_gibc_base(var.var_gibc)
  {
    InitDefaultValues();
  }
  
  
  inline void VarAxisymProblem::SetWaveVector(const R3& k)
  {
    kwave3D = k;
  }
  
  
  //! updates the wave vector with the new pulsation
  inline void VarAxisymProblem::UpdateWaveVector(const Real_wp& omega)
  {
    Real_wp coef = omega / Norm2(kwave3D);
    kwave3D *= coef;
  }
  

  //! retrieving 3-D wave vector
  inline const R3& VarAxisymProblem::GetWaveVector() const
  {
    return kwave3D;
  }
  
  
  //! retrieving phase origin (point where the phase is equal to 0)
  inline const R3& VarAxisymProblem::GetPhaseOrigin() const
  {
    return origine_phase3D;
  }

  
  inline void VarAxisymProblem::SetPhaseOrigin(const R3& pt)
  {
    origine_phase3D = pt;
  }

  
  inline const Complex_wp& VarAxisymProblem::GetWaveNumberInfinity() const
  {
    return kwave_envelope;
  }
  
  
  inline int VarAxisymProblem::GetNbModesSource() const
  {
    return 1;
  }


  //! returns exp(-i m theta) for complex numbers
  inline void VarAxisymProblem
  ::GetFourierMode(const Real_wp& teta, Complex_wp& val) const
  {
    val = exp(-Iwp*Real_wp(var_boundary.GetCurrentModeNumber())*teta);
  }


  //! returns cos(m theta) for real numbers
  inline void VarAxisymProblem
  ::GetFourierMode(const Real_wp& teta, Real_wp& val) const
  { 
    val = cos(var_boundary.GetCurrentModeNumber()*teta); 
  }

  
  inline bool VarAxisymProblem::IsVertexOnAxis(int n) const
  {
    return Vertex_On_Axe(n);
  }

  
  //! returns true if the element iquad is close to the axis r = 0
  inline bool VarAxisymProblem::IsElementNearAxis(int iquad) const
  {
    return Element_On_Axe(iquad);
  }
  
  
  inline bool VarAxisymProblem::NumberOfModesToBeComputed() const
  {
    return number_mode_to_be_computed;
  }
  
  
  inline Real_wp VarAxisymProblem::GetModeThreshold() const
  {
    return threshold_mode;
  }

  
  //! returns value of bessel function if pre-computed
  inline Real_wp VarAxisymProblem
  ::GetBessel_Value(int n, int num_point) const
  {
    if (n >= 0)
      return Jn_precomputed(n)(num_point);
    else
      {
	if (-n%2 == 0)
	  return Jn_precomputed(-n)(num_point);
	else
	  return -Jn_precomputed(-n)(num_point);
      }
  }


  template<class TypeEquation>
  inline VarComputationRCS_Axi::VarComputationRCS_Axi(EllipticProblem<TypeEquation>& var)
    : VarComputationRCS_Base<Dimension2>(var), var_problem(var), var_axisym(var)
  {
  }
  
  

#endif
    
}

#define MONTJOIE_FILE_VAR_AXISYM_PROBLEM_INLINE_CXX
#endif
