#ifndef MONTJOIE_FILE_DEFINE_SOURCE_MAXWELL_AXI_INLINE_CXX

namespace Montjoie
{

  /******************************
   * VirtualMaxwellAxiSourceFEM *
   ******************************/
  
  
  template<class T> template<class TypeEquation>
  inline VirtualMaxwellAxiSourceFEM<T>::VirtualMaxwellAxiSourceFEM(const EllipticProblem<TypeEquation>& var)
    : VirtualSourceFEM<T, Dimension2>(var)
  {
  }
  
  
  template<class T>
  inline bool VirtualMaxwellAxiSourceFEM<T>::IsNonNull_SourceH(const VectR2& s)
  {
    return false;
  }
  

  template<class T>
  inline void VirtualMaxwellAxiSourceFEM<T>::EvaluateSource_H(int i, int j, const R2& x, Vector<T>& f)
  {
  }

  
  /***********************
   * PlaneWaveMaxwellAxi *
   ***********************/
  

  template<class TypeEquation>
  inline PlaneWaveMaxwellAxi::PlaneWaveMaxwellAxi(const EllipticProblem<TypeEquation>& var)
    : var_maxwell(var)
  {
    InitDefaultValues();
  }
    

  template<class TypeEquation>
  inline IncidentWaveProjector_MaxwellAxi
  ::IncidentWaveProjector_MaxwellAxi(const EllipticProblem<TypeEquation>& var,
				     IncidentWaveField<Complex_wp, Dimension2>& u_inc)
    : IncidentWaveProjector<Complex_wp, Dimension2>(var, u_inc), PlaneWaveMaxwellAxi(var)
  {
  }
  
  
  /***********************************
   * DiffractedWaveSource_MaxwellAxi *
   ***********************************/
  
  
  //! constructor with given problem
  /*!
    \param[in] var given problem
   */
  template<class TypeEquation>
  inline DiffractedWaveSource_MaxwellAxi
  ::DiffractedWaveSource_MaxwellAxi(const EllipticProblem<TypeEquation>& var,
				    IncidentWaveField<Complex_wp, Dimension2>&)
    : VirtualMaxwellAxiSourceFEM<Complex_wp>(var), PlaneWaveMaxwellAxi(var), var_maxwell(var)
  {
    ref_domain = 1;
    pml_element = false;
  }
  
  
  template<class TypeEquation>
  inline DiffractedWaveSource<HarmonicMaxwellEquation_HcurlAxi>
  ::DiffractedWaveSource(const EllipticProblem<TypeEquation>& var,
			 IncidentWaveField<Complex_wp, Dimension2>& u_inc)
    : DiffractedWaveSource_MaxwellAxi(var, u_inc)
  {
  }

  
  /******************************
   * TotalWaveSource_MaxwellAxi *
   ******************************/


  //! constructor with given problem
  /*!
    \param[in] var given problem
   */
  template<class TypeEquation>
  inline TotalWaveSource_MaxwellAxi
  ::TotalWaveSource_MaxwellAxi(const EllipticProblem<TypeEquation>& var,
			       IncidentWaveField<Complex_wp, Dimension2>&)
    : VirtualMaxwellAxiSourceFEM<Complex_wp>(var), PlaneWaveMaxwellAxi(var), var_maxwell(var)
  {
    ref_domain = 1;
    pml_element = false;
  }


  template<class TypeEquation>
  inline TotalWaveSource<HarmonicMaxwellEquation_HcurlAxi>
  ::TotalWaveSource(const EllipticProblem<TypeEquation>& var,
		    IncidentWaveField<Complex_wp, Dimension2>& u_inc)
    : TotalWaveSource_MaxwellAxi(var, u_inc)
  {
  }


  /*******************************
   * VolumetricSource_MaxwellAxi *
   *******************************/
  
  
  //! default constructor with given problem
  template<class T> template<class TypeEquation>
  inline VolumetricSource_MaxwellAxi<T>
  ::VolumetricSource_MaxwellAxi(const EllipticProblem<TypeEquation>& var,
				const Vector<VectString>& param) 
    : VirtualMaxwellAxiSourceFEM<T>(var), var_problem(var), var_source(var),
      var_axi(var), var_boundary(var)
  { 
    InitDefaultValues(var.modified_formulation);
  }
  

  template<class TypeEquation>
  inline VolumetricSource<HarmonicMaxwellEquationAxiDG>
  ::VolumetricSource(const EllipticProblem<TypeEquation>& var,
		     const Vector<VectString>& param)
    : VolumetricSource_MaxwellAxi<Complex_wp>(var, param)
  {
  }


  template<class TypeEquation>
  inline VolumetricSource<HarmonicMaxwellEquation_HcurlAxi>
  ::VolumetricSource(const EllipticProblem<TypeEquation>& var,
		     const Vector<VectString>& param)
    : VolumetricSource_MaxwellAxi<Complex_wp>(var, param)
  {
  }  
}

#define MONTJOIE_FILE_DEFINE_SOURCE_MAXWELL_AXI_INLINE_CXX
#endif
