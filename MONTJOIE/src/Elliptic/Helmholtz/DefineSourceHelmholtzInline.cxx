#ifndef MONTJOIE_FILE_DEFINE_SOURCE_HELMHOLTZ_INLINE_CXX

namespace Montjoie
{

  //! constructor with given problem
  template<class Dimension> template<class TypeEquation>
  inline DiffractedWaveSource_Helm<Dimension>
  ::DiffractedWaveSource_Helm(const EllipticProblem<TypeEquation>& var,
			      IncidentWaveField<Complex_wp, Dimension>& u_inc)
    : VirtualSourceFEM<Complex_wp, Dimension>(var), incident_wave(u_inc),
      var_problem(var), var_helm(var), var_boundary(var),
      fct_imped(var.GetGenericImpedanceFunction())
  {
    InitDefaultValues();
  }


  //! constructor with a given problem
  template<class Dimension> template<class TypeEquation>
  inline DiffractedWaveSource<HelmholtzEquation<Dimension> >
  ::DiffractedWaveSource(const EllipticProblem<TypeEquation>& vars_helm,
                         IncidentWaveField<Complex_wp, Dimension>& f_inc)
    : DiffractedWaveSource_Helm<Dimension>(vars_helm, f_inc)
  {
  }
  
  
  //! constructor with given problem
  template<class Dimension> template<class TypeEquation>
  inline TotalWaveSource_Helm<Dimension>
  ::TotalWaveSource_Helm(const EllipticProblem<TypeEquation>& var,
			  IncidentWaveField<Complex_wp, Dimension>& u_inc)
    : VirtualSourceFEM<Complex_wp, Dimension>(var), incident_wave(u_inc),
      var_problem(var), var_helm(var), var_boundary(var),
      fct_imped_abc(var.GetAbsorbingImpedanceFunction())
  {
    this->coef = var_helm.mu0;
    this->omega = var_problem.GetOmega();
    this->element_insidePML = false;
    this->rho0 = var_helm.rho0; this->mu0 = var_helm.mu0;
  }
  

  //! constructor with a given problem
  template<class Dimension> template<class TypeEquation>
  inline TotalWaveSource<HelmholtzEquation<Dimension> >
  ::TotalWaveSource(const EllipticProblem<TypeEquation>& vars_helm,
		    IncidentWaveField<Complex_wp, Dimension>& f_inc)
    : TotalWaveSource_Helm<Dimension>(vars_helm, f_inc) 
  {
  }


  //! constructor with a given problem
  template<class Dimension> template<class TypeEquation>
  inline DiffractedWaveSource<HelmholtzEquationDG<Dimension> >
  ::DiffractedWaveSource(const EllipticProblem<TypeEquation>& vars_helm,
			 IncidentWaveField<Complex_wp, Dimension>& f_inc)
    : DiffractedWaveSource_HelmDG<Dimension>(vars_helm, f_inc) 
  {
  }


  //! constructor with a given problem
  template<class Dimension> template<class TypeEquation>
  inline TotalWaveSource<HelmholtzEquationDG<Dimension> >
  ::TotalWaveSource(const EllipticProblem<TypeEquation>& vars_helm,
		    IncidentWaveField<Complex_wp, Dimension>& f_inc)
    : TotalWaveSource_HelmDG<Dimension>(vars_helm, f_inc)
  {
  }


  //! constructor with given problem
  template<class Dimension> template<class TypeEquation>
  inline DiffractedWaveSource_HelmDG<Dimension>
  ::DiffractedWaveSource_HelmDG(const EllipticProblem<TypeEquation>& var,
				IncidentWaveField<Complex_wp, Dimension>& u_inc)
    : VirtualSourceFEM<Complex_wp, Dimension>(var), incident_wave(u_inc),
      var_problem(var), var_helm(var), var_boundary(var)
  {
    InitDefaultValues();
  }


  //! constructor with given problem
  template<class Dimension> template<class TypeEquation>
  inline TotalWaveSource_HelmDG<Dimension>
  ::TotalWaveSource_HelmDG(const EllipticProblem<TypeEquation>& var,
			   IncidentWaveField<Complex_wp, Dimension>& u_inc)
    : VirtualSourceFEM<Complex_wp, Dimension>(var), incident_wave(u_inc),
      var_problem(var), var_helm(var), var_boundary(var)
  {
    InitDefaultValues();
  }


  template<class Dimension> template<class TypeEquation>
  inline IncidentWaveProjector_HelmDG<Dimension>
  ::IncidentWaveProjector_HelmDG(const EllipticProblem<TypeEquation>& var,
				 IncidentWaveField<Complex_wp, Dimension>& u_inc)
    : IncidentWaveProjector<Complex_wp, Dimension>(var, u_inc), var_helm(var), var_problem(var)
  {
  }


  template<class T, class Dimension> template<class TypeEquation>
  inline VolumetricSource_HelmDG<T, Dimension>
  ::VolumetricSource_HelmDG(const EllipticProblem<TypeEquation>& var,
			    const Vector<VectString>& param)
    : VolumetricSource_Base<T, Dimension>(var, param), var_helm(var), var_problem(var), var_boundary(var)
  {
    
  }

  template<class Dimension>
  inline VolumetricSource<HelmholtzEquationDG<Dimension> >
  ::VolumetricSource(const EllipticProblem<HelmholtzEquationDG<Dimension> >& var,
		     const Vector<VectString>& param)
    : VolumetricSource_HelmDG<Complex_wp, Dimension>(var, param)
  {
  }


}

#define MONTJOIE_FILE_DEFINE_SOURCE_HELMHOLTZ_INLINE_CXX
#endif
