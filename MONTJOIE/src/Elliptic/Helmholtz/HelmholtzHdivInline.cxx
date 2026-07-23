#ifndef MONTJOIE_FILE_HELMHOLTZ_HDIV_INLINE_CXX

namespace Montjoie
{

  template<class Dimension> template<class TypeEquation>
  inline VarHelmholtzHdiv_Base<Dimension>::VarHelmholtzHdiv_Base(EllipticProblem<TypeEquation>& var)
    : var_problem(var)
  {
  }
  

  //! returns the number of different physical media
  template<class Dimension>
  inline int VarHelmholtzHdiv_Base<Dimension>::GetNbPhysicalIndices() const
  {
    return ref_rho0.GetM();
  }
  
  
  //! Are the physical indices variable inside element i
  template<class Dimension>
  inline bool VarHelmholtzHdiv_Base<Dimension>::IsVaryingMedia(int i) const
  {
    return true;
  }


  //! returns the velocity of waves in a media
  template<class Dimension>
  inline Real_wp VarHelmholtzHdiv_Base<Dimension>::GetVelocityOfMedia(int ref) const
  {
    // we do not take into account the flow
    return this->ref_c0(ref).GetConstant();
  }
  

  //! returns the velocity of waves at infinity
  template<class Dimension>
  inline Real_wp VarHelmholtzHdiv_Base<Dimension>::GetVelocityOfInfinity() const
  {
    return 1.0;
  }


  //! constructor with given problem
  template<class Dimension> template<class TypeEquation>
  inline DiffractedWaveSource_HelmHdiv<Dimension>::
  DiffractedWaveSource_HelmHdiv(const EllipticProblem<TypeEquation>& vars_helm,
				IncidentWaveField<Complex_wp, Dimension>& u_inc)
    : VirtualSourceFEM<Complex_wp, Dimension>(vars_helm), incident_wave(u_inc), var_problem(vars_helm)
  {
    InitDefaultValues();
  }


  template<class Dimension> template<class TypeEquation>
  inline IncidentWaveProjector_HelmHdiv<Dimension>
  ::IncidentWaveProjector_HelmHdiv(const EllipticProblem<TypeEquation>& var,
				   IncidentWaveField<Complex_wp, Dimension>& u_inc)
    : IncidentWaveProjector<Complex_wp, Dimension>(var, u_inc), var_helm(var), var_problem(var)
  {    
  }  


  //! constructor with a given problem
  template<class Dimension> template<class TypeEquation>
  inline DiffractedWaveSource<HelmholtzEquationHdiv<Dimension> >
  ::DiffractedWaveSource(const EllipticProblem<TypeEquation>& vars_helm,
			      IncidentWaveField<Complex_wp, Dimension>& u_inc)
    : DiffractedWaveSource_HelmHdiv<Dimension>(vars_helm, u_inc)
  {
  }
  

  template<class Dimension> template<class TypeEquation>
  inline ImpedanceABC<Complex_wp, HelmholtzEquationHdiv<Dimension> >
  ::ImpedanceABC(const EllipticProblem<TypeEquation>& var)
    : ImpedanceFunction_Base<Complex_wp, Dimension>(var), var_helm(var), var_problem(var)
  {
    order = var.GetOrderAbsorbingCondition();
    delta = var.gamma_cla_coef;
  }
  


}

#define MONTJOIE_FILE_HELMHOLTZ_HDIV_INLINE_CXX
#endif
