#ifndef MONTJOIE_FILE_VAR_HARMONIC_BASE_HXX

namespace Montjoie
{
  
  //! base class for time-harmonic (or stationary) problems
  /*!
    This class regroups generic methods used to solve a stationary problem
    with finite elements, mainly computation of the matrix, computation of right hand sides
    and outputs of the solution.
  */
  template<class Complexe, class Dimension>
  class VarHarmonic_Base : public VarProblem<Dimension>,
			   public VarBoundaryCondition<Complexe, Dimension>,
			   public VarOutputProblem<Dimension>,
			   public VarSourceProblem<Complexe, Dimension>
  {    

  private:
    VarGeneralizedImpedance_Base& var_gibc_base;

#ifdef MONTJOIE_WITH_TRANSMISSION
    VarTransmission_Base<Dimension>& var_transmission_base;
#endif

    VarComputationRCS_Base<Dimension>& output_rcs_base;

    typedef Vector<Complexe> VectComplexe;
    typedef typename Dimension::VectR_N VectR_N;
    
  public :

    template<class TypeEquation>
    VarHarmonic_Base(EllipticProblem<TypeEquation>& var);
    
    // Inline methods
    virtual bool IsComplexProblem() const;
    
    // other methods
    virtual void SetInputData(const string& description_field, const VectString& parameters);    
    
    virtual void PerformOtherInitializations();    
    virtual void Restart();
    
    // method to run global simulation
    void ConstructAll(const string& input_file, const string& name_elt, const string&,
		      All_LinearSolver*&, bool compute_rho = true,
		      bool delete_points = true, int num = -1);

    void ConstructAll(const Vector<string>& lines_data_file, const string& input_file,
		      const string&, const string&, All_LinearSolver*&,
		      bool compute_rho = true, bool delete_points = true, int num = -1);
    
    virtual void RunAll(const string& input_file, const string& name_elt, const string&, int num = -1);
    virtual void RunAll(const Vector<string>& lines_data_file, const string& input_file,
			const string&, const string&, int num = -1);
    
    // other methods
    virtual void GetMemoryUsed(map<string, size_t>& var) const;
    virtual size_t GetMemorySize() const;

    void GetInverseSquareRootMassMatrix(VectComplexe& Dh);
    virtual void GetMassMatrix(VectComplexe& Dh, bool assemble = true) const = 0;

    void CopyInputData(const VarProblem_Base& var_base);
    
    // Inline methods
    void ComputeMassMatrix(bool compute_rho = true, bool delete_points = true);    
    void ComputeStoreCoefficientsPML(int i1, int ne, const VectR_N& AllPoints);

    void SetComputationFarPoints(VectR_N& Points, Real_wp dt);

    virtual TransparencySolver_Base* GetNewTransparentSolver(All_LinearSolver& glob_solver) = 0;    
    virtual EigenProblemMontjoie<Complexe>* GetNewEigenSolver(All_LinearSolver& glob_solver) = 0;
    virtual EigenProblemMontjoie<Real_wp>* GetNewEigenSolver(All_LinearSolver& glob_solver, const Real_wp& x) = 0;
    virtual EigenProblemMontjoie<Complex_wp>* GetNewEigenSolver(All_LinearSolver& glob_solver, const Complex_wp& x) = 0;
    virtual PolynomialEigenProblemMontjoie<Complexe>* GetNewPolynomialEigenSolver(All_LinearSolver& glob_solver) = 0;
    
    void InitVarGrid();
    void FindElementsInsidePML();

    void TreatGibc(const IVect& Epart);
    void InitGibcReferences(int N);

#ifdef MONTJOIE_WITH_TRANSMISSION
    void PartMeshTransmission();
    void TreatTransmission(const IVect& epart);
    void SendTransmissionDofs(const IVect& num, int& nb0, int& nb1, IVect& info0, IVect& info1);
    void DistributeTransmissionDofs(int nb0, int nb1, IVect& info0, IVect& info1);
#endif
    
  };
} 

#define MONTJOIE_FILE_VAR_HARMONIC_BASE_HXX
#endif
