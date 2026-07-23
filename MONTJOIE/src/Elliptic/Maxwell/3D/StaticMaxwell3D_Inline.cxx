#ifndef MONTJOIE_FILE_STATIC_MAXWELL_3D_INLINE_CXX

namespace Montjoie
{

  template<class TypeEquation>
  inline IncidentWaveProjector_StaticMaxwell3D
  ::IncidentWaveProjector_StaticMaxwell3D(const EllipticProblem<TypeEquation>& var,
                                          IncidentWaveField<Real_wp, Dimension3>& u_inc)
    : IncidentWaveProjector<Real_wp, Dimension3>(var, u_inc), var_source(var),
      var_boundary(var), var_problem(var)     
  {
  }
    

  //! default constructor
  template<class TypeEquation>
  inline DiffractedWaveSource_StaticMaxwell3D
  ::DiffractedWaveSource_StaticMaxwell3D(const EllipticProblem<TypeEquation>& var,
                                         IncidentWaveField<Real_wp, Dimension3>& u_inc)
    : VirtualSourceFEM<Real_wp, Dimension3>(var), incident_wave(u_inc), var_source(var),
      var_boundary(var), var_problem(var), var_maxwell(var)
  {
    scalar_eq = true;
    dirichlet_cond = true;
  }    
  

  //! constructor with given problem
  template<class TypeEquation>
  inline DiffractedWaveSource<StaticMaxwellEquation_3D>
  ::DiffractedWaveSource(const EllipticProblem<TypeEquation>& vars_helm,
			 IncidentWaveField<Real_wp, Dimension3>& fsrc_)
    : DiffractedWaveSource_StaticMaxwell3D(vars_helm, fsrc_) 
  {
  }
  

  //! default constructor
  template<class TypeEquation>
  inline TotalWaveSource_StaticMaxwell3D
  ::TotalWaveSource_StaticMaxwell3D(const EllipticProblem<TypeEquation>& var,
                                    IncidentWaveField<Real_wp, Dimension3>& u_inc)
    : VirtualSourceFEM<Real_wp, Dimension3>(var), incident_wave(u_inc), var_source(var),
      var_boundary(var), var_problem(var), var_maxwell(var)
  {
    scalar_eq = true;
  }    
  

  //! constructor with given problem
  template<class TypeEquation>
  inline TotalWaveSource<StaticMaxwellEquation_3D>
  ::TotalWaveSource(const EllipticProblem<TypeEquation>& vars_helm,
		    IncidentWaveField<Real_wp, Dimension3>& fsrc)
    : TotalWaveSource_StaticMaxwell3D(vars_helm, fsrc) 
  {
  }

  template<class Complexe> template<class TypeEquation>
  inline ImpedanceABC<Complexe, StaticMaxwellEquation_3D>::ImpedanceABC(const EllipticProblem<TypeEquation>& var)
    : ImpedanceABC_Maxwell3D<Complexe, Real_wp>(var)
  {
  }
  
  template<class T>
  inline CondensationBlockSolver<T, StaticMaxwellEquation_3D >
  ::CondensationBlockSolver(EllipticProblem<StaticMaxwellEquation_3D >& var)
    : CondensationBlockSolver_Maxwell3D<T, StaticMaxwellEquation_3D >(var)
  {
  }

  template<class T>
  inline CondensationBlockSolver<T, StaticMaxwellEquationHdg_3D>
  ::CondensationBlockSolver(EllipticProblem<StaticMaxwellEquationHdg_3D>& var)
    : CondensationBlockSolver_MaxwellHdg3D<T, Real_wp>(var)
  {
  }

  inline FemMatrixFreeClass<Real_wp, StaticMaxwellEquation_3D>
  ::FemMatrixFreeClass(const EllipticProblem<StaticMaxwellEquation_3D>& var_)
    : FemMatrixFreeClass_Eq<Real_wp, StaticMaxwellEquation_3D>(var_)
  {
  }

  
  inline FemMatrixFreeClass<Complex_wp, StaticMaxwellEquation_3D>
  ::FemMatrixFreeClass(const EllipticProblem<StaticMaxwellEquation_3D>& var_)
    : FemMatrixFreeClass_Eq<Complex_wp, StaticMaxwellEquation_3D>(var_)
  {
  }

  inline VolumetricSource<StaticMaxwellEquationHdg_3D>
  ::VolumetricSource(const EllipticProblem<StaticMaxwellEquationHdg_3D>& var,
		     const Vector<VectString>& param)
    : VolumetricSource_MaxwellHdg3D(var, param)
  {
  }
  
  inline FemMatrixFreeClass<Real_wp, StaticMaxwellEquationHdg_3D>
  ::FemMatrixFreeClass(const EllipticProblem<StaticMaxwellEquationHdg_3D>& var_)
    : FemMatrixFreeClass_Eq<Real_wp, StaticMaxwellEquationHdg_3D>(var_)
  {
  }

  
  inline FemMatrixFreeClass<Complex_wp, StaticMaxwellEquationHdg_3D>
  ::FemMatrixFreeClass(const EllipticProblem<StaticMaxwellEquationHdg_3D>& var_)
    : FemMatrixFreeClass_Eq<Complex_wp, StaticMaxwellEquationHdg_3D>(var_)
  {
  }
}

#define MONTJOIE_FILE_STATIC_MAXWELL_3D_INLINE_CXX
#endif
