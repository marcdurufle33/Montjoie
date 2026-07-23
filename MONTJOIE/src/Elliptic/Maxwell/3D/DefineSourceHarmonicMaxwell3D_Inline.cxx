#ifndef MONTJOIE_FILE_DEFINE_SOURCE_HARMONIC_MAXWELL_3D_INLINE_CXX

namespace Montjoie
{

  //! constructor with given problem
  template<class TypeEquation>
  inline DiffractedWaveSource_Maxwell3D
  ::DiffractedWaveSource_Maxwell3D(const EllipticProblem<TypeEquation>& var,
				   IncidentWaveField<Complex_wp, Dimension3>& u_inc)
    : VirtualSourceFEM<Complex_wp, Dimension3>(var), incident_wave(u_inc),
      var_problem(var), var_maxwell(var),
      fct_imped(var.GetGenericImpedanceFunction())
  {
    InitDefaultValues();
  }
  
  
  //! constructor with a given problem
  template<class TypeEquation>
  inline DiffractedWaveSource<HarmonicMaxwellEquation_3D >
  ::DiffractedWaveSource(const EllipticProblem<TypeEquation>& vars_helm,
                         IncidentWaveField<Complex_wp, Dimension3>& f_inc)
    : DiffractedWaveSource_Maxwell3D(vars_helm, f_inc)
  {
  }
  
  
  //! constructor with given problem
  template<class TypeEquation>
  inline TotalWaveSource_Maxwell3D
  ::TotalWaveSource_Maxwell3D(const EllipticProblem<TypeEquation>& var,
			      IncidentWaveField<Complex_wp, Dimension3>& u_inc)
    : VirtualSourceFEM<Complex_wp, Dimension3>(var), incident_wave(u_inc),
      var_problem(var), var_maxwell(var), var_boundary(var)
      //, fct_imped_abc(var.GetAbsorbingImpedanceFunction())
  {
    InitDefaultValues();
  }
  

  //! constructor with a given problem
  template<class TypeEquation>
  inline TotalWaveSource<HarmonicMaxwellEquation_3D >
  ::TotalWaveSource(const EllipticProblem<TypeEquation>& vars_helm,
		    IncidentWaveField<Complex_wp, Dimension3>& f_inc)
    : TotalWaveSource_Maxwell3D(vars_helm, f_inc) 
  {
  }


  //! constructor with a given problem
  template<class TypeEquation>
  inline DiffractedWaveSource<HarmonicMaxwellEquation_3D_DG>
  ::DiffractedWaveSource(const EllipticProblem<TypeEquation>& var,
			 IncidentWaveField<Complex_wp, Dimension3>& u_inc)
    : VirtualSourceFEM<Complex_wp, Dimension3>(var), incident_wave(u_inc),
      var_problem(var), var_maxwell(var)
      //fct_imped(var.GetGenericImpedanceFunction())
  {
    InitDefaultValues();
  }


  //! constructor with a given problem
  template<class TypeEquation>
  inline TotalWaveSource<HarmonicMaxwellEquation_3D_DG>
  ::TotalWaveSource(const EllipticProblem<TypeEquation>& var,
		    IncidentWaveField<Complex_wp, Dimension3>& u_inc)
    : VirtualSourceFEM<Complex_wp, Dimension3>(var), incident_wave(u_inc),
      var_problem(var)
      //fct_imped_abc(var.GetAbsorbingImpedanceFunction())
  {
    InitDefaultValues();
  }
  

  template<class TypeEquation>
  inline IncidentWaveProjector_Maxwell3D_DG
  ::IncidentWaveProjector_Maxwell3D_DG(const EllipticProblem<TypeEquation>& var,
				       IncidentWaveField<Complex_wp, Dimension3>& u_inc)
    : IncidentWaveProjector<Complex_wp, Dimension3>(var, u_inc), var_problem(var)
  {
    var_problem.GetPolarization(E_0);    
    TimesProd(var_problem.GetWaveVector(), E_0, H_0);
  }


  template<class T> template<class TypeEquation>
  inline VolumetricSource_MaxwellHdg3D<T>
  ::VolumetricSource_MaxwellHdg3D(const EllipticProblem<TypeEquation>& var,
                                const Vector<VectString>& param)
    : VolumetricSource_Base<T, Dimension3>(var, param),
      var_maxwell(var), var_problem(var), var_boundary(var)
  {
  }


  inline VolumetricSource<HarmonicMaxwellEquationHdg_3D>
  ::VolumetricSource(const EllipticProblem<HarmonicMaxwellEquationHdg_3D>& var,
		     const Vector<VectString>& param)
    : VolumetricSource_MaxwellHdg3D<Complex_wp>(var, param)
  {
  }
  
}

#define MONTJOIE_FILE_DEFINE_SOURCE_HARMONIC_MAXWELL_3D_INLINE_CXX
#endif
