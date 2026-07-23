#ifndef MONTJOIE_FILE_VAR_HARMONIC_INLINE_CXX

namespace Montjoie
{
    
  //! returns leaf class
  template<class TypeEquation>
  inline EllipticProblem<TypeEquation>& VarHarmonic<TypeEquation>::GetLeafClass()
  {
    return static_cast<EllipticProblem<TypeEquation>& >(*this);
  }

  
  //! returns leaf class
  template<class TypeEquation>
  inline const EllipticProblem<TypeEquation>& VarHarmonic<TypeEquation>::GetLeafClass() const
  {
    return static_cast<const EllipticProblem<TypeEquation>& >(*this);
  }

  
  //! allocation of projector for the incident field
  template<class TypeEquation>
  inline IncidentWaveProjector<typename TypeEquation::Complexe, typename TypeEquation::Dimension>*
  VarHarmonic<TypeEquation>
  ::GetNewIncidentProjector(int n, const Vector<VectString>& param,
			    IncidentWaveField<Complexe, Dimension>& u_inc) const
  {
    return new IncidentWaveProjector<Complexe, Dimension>(GetLeafClass(), u_inc);
  }

  
  //! allocation of a modal source
  template<class TypeEquation>
  inline ModalSourceBoundary_Dim<typename TypeEquation::Complexe, typename TypeEquation::Dimension>* 
  VarHarmonic<TypeEquation>::GetNewModalSourceEquation() const
  {
    return new ModalSourceBoundary_Dim<Complexe, Dimension>(GetLeafClass());
  }
  

  template<class TypeEquation>
  inline ImpedanceGeneric<typename TypeEquation::Complexe, TypeEquation>& VarHarmonic<TypeEquation>
  ::GetGenericImpedanceFunction() const
  {
    ImpedanceGeneric<Complexe, TypeEquation>& fct
      = const_cast<ImpedanceGeneric<Complexe, TypeEquation>& >(fct_impedance_generic);
    
    return fct;
  }
  
  
  template<class TypeEquation>
  inline ImpedanceABC<typename TypeEquation::Complexe, TypeEquation>& VarHarmonic<TypeEquation>
  ::GetAbsorbingImpedanceFunction() const
  {
    ImpedanceABC<Complexe, TypeEquation>& fct
      = const_cast<ImpedanceABC<Complexe, TypeEquation>& >(fct_impedance_absorbing);
    
    return fct;
  }

  
  template<class TypeEquation>
  inline ImpedanceFunction_Base<Real_wp, typename TypeEquation::Dimension>* 
  VarHarmonic<TypeEquation>::GetNewImpedanceABC(const Real_wp&) const
  {
    return new ImpedanceABC<Real_wp, TypeEquation>(GetLeafClass());
  }
  
  
  template<class TypeEquation>
  inline ImpedanceFunction_Base<Complex_wp, typename TypeEquation::Dimension>* 
  VarHarmonic<TypeEquation>::GetNewImpedanceABC(const Complex_wp&) const
  {
    return new ImpedanceABC<Complex_wp, TypeEquation>(GetLeafClass()); 
  }
  
  
  template<class TypeEquation>
  inline ImpedanceFunction_Base<Real_wp, typename TypeEquation::Dimension>* 
  VarHarmonic<TypeEquation>::GetNewImpedanceGeneric(const Real_wp&) const
  { 
    return new ImpedanceGeneric<Real_wp, TypeEquation>(GetLeafClass());
  }
  
  
  template<class TypeEquation>
  inline ImpedanceFunction_Base<Complex_wp, typename TypeEquation::Dimension>* 
  VarHarmonic<TypeEquation>::GetNewImpedanceGeneric(const Complex_wp&) const
  { 
    return new ImpedanceGeneric<Complex_wp, TypeEquation>(GetLeafClass()); 
  }
  

#ifdef MONTJOIE_WITH_HIGH_CONDUCTIVITY_MODEL    
  template<class TypeEquation>
  inline ImpedanceFunction_Base<Real_wp, typename TypeEquation::Dimension>* 
  VarHarmonic<TypeEquation>::GetNewImpedanceHighConductivity(const Real_wp&) const
  { 
    return new ImpedanceHighConductivity<Real_wp, TypeEquation>(GetLeafClass());
  }

  
  template<class TypeEquation>
  inline ImpedanceFunction_Base<Complex_wp, typename TypeEquation::Dimension>* 
  VarHarmonic<TypeEquation>::GetNewImpedanceHighConductivity(const Complex_wp&) const
  { 
    return new ImpedanceHighConductivity<Complex_wp, TypeEquation>(GetLeafClass());
  }
#endif

  
  template<class TypeEquation>
  inline FemMatrixFreeClass_Base<Real_wp>* VarHarmonic<TypeEquation>
  ::GetNewIterativeMatrix(const Real_wp&) const
  {
    FemMatrixFreeClass_Base<Real_wp>* Mh = new FemMatrixFreeClass<Real_wp, TypeEquation>(GetLeafClass());
    Mh->Reallocate(this->GetNbDof(), this->GetNbDof());
    return Mh;
  }
  
  
  template<class TypeEquation>
  inline FemMatrixFreeClass_Base<Complex_wp>* VarHarmonic<TypeEquation>
  ::GetNewIterativeMatrix(const Complex_wp&) const
  {
    FemMatrixFreeClass_Base<Complex_wp>* Mh = new FemMatrixFreeClass<Complex_wp, TypeEquation>(GetLeafClass());
    Mh->Reallocate(this->GetNbDof(), this->GetNbDof());
    return Mh;
  }
  
  
  template<class TypeEquation>
  inline DistributedProblem<typename TypeEquation::Dimension>* VarHarmonic<TypeEquation>
  ::GetNewEllipticProblem() const
  {
    return new EllipticProblem<TypeEquation>();
  }
  
  
  template<class TypeEquation>
  inline All_LinearSolver* VarHarmonic<TypeEquation>::GetNewLinearSolver()
  {
    return new All_LinearSolver(GetLeafClass());
  }

  
  //! returns a new preconditioning class
  template<class TypeEquation>
  inline All_Preconditioner_Base<Real_wp>* 
  VarHarmonic<TypeEquation>::GetNewPreconditioning(const Real_wp&)
  {
    return new All_Preconditioner<Real_wp, TypeEquation>(GetLeafClass());
  }
  
  
  //! returns a new preconditioning class
  template<class TypeEquation>
  inline All_Preconditioner_Base<Complex_wp>*
  VarHarmonic<TypeEquation>::GetNewPreconditioning(const Complex_wp&)
  {
    return new All_Preconditioner<Complex_wp, TypeEquation>(GetLeafClass());
  }
  
  
  template<class TypeEquation>
  inline EigenProblemMontjoie<typename TypeEquation::Complexe>* VarHarmonic<TypeEquation>
  ::GetNewEigenSolver(All_LinearSolver& glob_solver)
  {
    return new EigenProblemMontjoie<Complexe>(GetLeafClass(), glob_solver);
  }
  
  
  template<class TypeEquation>
  inline EigenProblemMontjoie<Real_wp>* VarHarmonic<TypeEquation>
  ::GetNewEigenSolver(All_LinearSolver& glob_solver, const Real_wp&)
  {
    return new EigenProblemMontjoie<Real_wp>(GetLeafClass(), glob_solver);
  }
  
  
  template<class TypeEquation>
  inline EigenProblemMontjoie<Complex_wp>* VarHarmonic<TypeEquation>
  ::GetNewEigenSolver(All_LinearSolver& glob_solver, const Complex_wp&)
  {
    return new EigenProblemMontjoie<Complex_wp>(GetLeafClass(), glob_solver);
  }


  template<class TypeEquation>
  inline PolynomialEigenProblemMontjoie<typename TypeEquation::Complexe>* VarHarmonic<TypeEquation>
  ::GetNewPolynomialEigenSolver(All_LinearSolver& glob_solver)
  {
    return new PolynomialEigenProblemMontjoie<Complexe>(GetLeafClass(), glob_solver);
  }

  template<class TypeEquation>
  inline TransparencySolver_Base* VarHarmonic<TypeEquation>
  ::GetNewTransparentSolver(All_LinearSolver& glob_solver)
  {
    return new TransparencySolver<TypeEquation>(GetLeafClass(), glob_solver);
  }

  
  template<class TypeEquation>
  inline CondensationBlockSolver_Base<Real_wp>* VarHarmonic<TypeEquation>
  ::GetNewCondensationSolver(const Real_wp&)
  {
    return new CondensationBlockSolver<Real_wp, TypeEquation>(GetLeafClass());
  }


  template<class TypeEquation>
  inline CondensationBlockSolver_Base<Complex_wp>* VarHarmonic<TypeEquation>
  ::GetNewCondensationSolver(const Complex_wp&)
  {
    return new CondensationBlockSolver<Complex_wp, TypeEquation>(GetLeafClass());
  }
  
}

#define MONTJOIE_FILE_VAR_HARMONIC_INLINE_CXX
#endif
