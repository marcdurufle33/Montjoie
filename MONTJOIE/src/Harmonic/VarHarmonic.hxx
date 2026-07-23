#ifndef MONTJOIE_FILE_VAR_HARMONIC_HXX

namespace Montjoie
{
  
  //! base class for time-harmonic (or stationary) problems
  /*!
    This class regroups generic methods used to solve a stationary problem
    with finite elements, mainly computation of the matrix, computation of right hand sides
    and outputs of the solution.
  */
  template<class TypeEquation>
  class VarHarmonic : public VarHarmonic_Base<typename TypeEquation::Complexe,
					      typename TypeEquation::Dimension>
  {
    
  public :
    //! the equation contains the dimension
    typedef typename TypeEquation::Dimension Dimension;    
    typedef typename Dimension::VectR_N VectR_N;
    
    typedef typename TypeEquation::Complexe Complexe; //!< complex or real
    typedef Vector<Complexe> VectComplexe; //!< "complex" vector
    
    ImpedanceABC<Complexe, TypeEquation> fct_impedance_absorbing;
    ImpedanceGeneric<Complexe, TypeEquation> fct_impedance_generic;
    
#ifdef MONTJOIE_WITH_HIGH_CONDUCTIVITY_MODEL    
    ImpedanceHighConductivity<Complexe, TypeEquation> fct_impedance_high_conduc;
#endif

    VarComputationRCS<TypeEquation> output_rcs_param;
    
#ifdef MONTJOIE_WITH_TRANSMISSION
    VarTransmission<TypeEquation> var_transmission;
#endif

    VarGeneralizedImpedance<TypeEquation> var_gibc;
    
  public :
    VarHarmonic();

    void GetMemoryUsed(map<string, size_t>& var) const;
    
    // returns a pointer to the source object
    VirtualSourceFEM<Complexe, Dimension>* GetNewSourceEquationObject(int n);
    
    // returns a pointer to a Dirac source object
    VirtualSourceFEM<Complexe, Dimension>* GetNewDiracSourceEquationObject(const Vector<VectString>& param) const;
    
    // mass matrix
    void GetMassMatrix(VectComplexe& Dh, bool assemble = true) const;    

    template<class TypeEq2>
    void CopyIndices(const EllipticProblem<TypeEq2>& var);

  protected:
    void ComputeArraySpecificEquation();

  public:
    // Inline methods
    EllipticProblem<TypeEquation>& GetLeafClass();
    const EllipticProblem<TypeEquation>& GetLeafClass() const;

    IncidentWaveProjector<Complexe, Dimension>*
    GetNewIncidentProjector(int n, const Vector<VectString>&,
			    IncidentWaveField<Complexe, Dimension>& u_inc) const;

    ModalSourceBoundary_Dim<Complexe, Dimension>* GetNewModalSourceEquation() const;
    
    ImpedanceGeneric<Complexe, TypeEquation>& GetGenericImpedanceFunction() const;
    ImpedanceABC<Complexe, TypeEquation>& GetAbsorbingImpedanceFunction() const;

    ImpedanceFunction_Base<Real_wp, Dimension>* 
    GetNewImpedanceABC(const Real_wp&) const;

    ImpedanceFunction_Base<Complex_wp, Dimension>* 
    GetNewImpedanceABC(const Complex_wp&) const;

    ImpedanceFunction_Base<Real_wp, Dimension>* 
    GetNewImpedanceGeneric(const Real_wp&) const;

    ImpedanceFunction_Base<Complex_wp, Dimension>* 
    GetNewImpedanceGeneric(const Complex_wp&) const;

#ifdef MONTJOIE_WITH_HIGH_CONDUCTIVITY_MODEL    
    ImpedanceFunction_Base<Real_wp, Dimension>* 
    GetNewImpedanceHighConductivity(const Real_wp&) const;

    ImpedanceFunction_Base<Complex_wp, Dimension>* 
    GetNewImpedanceHighConductivity(const Complex_wp&) const;
#endif

    FemMatrixFreeClass_Base<Real_wp>* GetNewIterativeMatrix(const Real_wp&) const;
    FemMatrixFreeClass_Base<Complex_wp>* GetNewIterativeMatrix(const Complex_wp&) const;

    DistributedProblem<Dimension>* GetNewEllipticProblem() const;

    All_LinearSolver* GetNewLinearSolver();

    All_Preconditioner_Base<Real_wp>* GetNewPreconditioning(const Real_wp&);
    All_Preconditioner_Base<Complex_wp>* GetNewPreconditioning(const Complex_wp&);
    
    EigenProblemMontjoie<Complexe>* GetNewEigenSolver(All_LinearSolver& glob_solver);
    EigenProblemMontjoie<Real_wp>* GetNewEigenSolver(All_LinearSolver& glob_solver, const Real_wp&);
    EigenProblemMontjoie<Complex_wp>* GetNewEigenSolver(All_LinearSolver& glob_solver, const Complex_wp&);
    PolynomialEigenProblemMontjoie<Complexe>* GetNewPolynomialEigenSolver(All_LinearSolver& glob_solver);
    
    TransparencySolver_Base* GetNewTransparentSolver(All_LinearSolver& glob_solver);

    virtual CondensationBlockSolver_Base<Real_wp>* GetNewCondensationSolver(const Real_wp&);
    virtual CondensationBlockSolver_Base<Complex_wp>* GetNewCondensationSolver(const Complex_wp&);
    
  };
  
}

#define MONTJOIE_FILE_VAR_HARMONIC_HXX  
#endif
