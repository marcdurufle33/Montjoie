#ifndef MONTJOIE_FILE_VAR_INSTATIONARY_HXX

#include "TimeSchemes.hxx"
 
namespace Montjoie
{

  //! parameters for mass matrix
  class MassMatrixParameter
  {
    friend class VarInstationary_Base;
    
  protected:
    bool store_mass_matrix;
    bool iterative_solver_mass;
    Real_wp stopping_criterion;
    int nb_max_iteration;
    int print_level;
    
  public:
    MassMatrixParameter();

    inline bool MassMatrixStored() const { return store_mass_matrix; }
    inline bool UseIterativeSolver() const { return iterative_solver_mass; }
    inline int GetNbMaxIterations() const { return nb_max_iteration; }
    inline Real_wp GetStoppingCriterion() const { return stopping_criterion; }
    inline int GetPrintLevel() const { return print_level; }
    
  };
    
  //! Base class for instationnary problems, time-domain maxwell equations, 
  //! wave equation, aeroacoustic equations ...
  class VarInstationary_Base : public VirtualOdeSystem<Real_wp>,
			       public InputDataProblem_Base
  {
  public :
    //! linear solver used for implicit time schemes
    All_LinearSolver* glob_solver;

    //! set of linear solvers
    Vector<All_LinearSolver*> vec_solver;
    
    //! stiffness matrix
    FemMatrixFreeClass_Base<Real_wp>* Glob_mat_Kh;

    //! mass matrix
    VirtualMassMatrix* Glob_mat_Dh;
    
    //! object interfacing the time schemes
    All_TimeScheme<Real_wp> var_time_scheme;
    
    //! stiffness matrices (if stored)
    DistributedMatrix<Real_wp, General, ArrayRowSparse> Glob_mat_RhS, Glob_mat_RhV, Glob_mat_Ch_Lambda;    
    
    bool load_reprise, save_reprise;
    int load_iter_reprise, save_iter_reprise, bascule_reprise;
    string path_reprise;
    //! parametres de reprise
    ParameterOutputReprise output_reprise;

  protected:
    bool time_step_to_be_computed; //!< true if the time step deltat has to be computed
    //!< number of neighboring elements to take into account for the cfl evaluation
    Real_wp deltat; //!< global time step
    Real_wp initial_time; //!< initial time t0
    Real_wp final_time; //!< final time tf
    Real_wp max_norme_sol;
    
    MatrixVectorProductLevel level_time_scheme;
    string file_with_local_dt;
    int global_level_operator;
    const int dim_N;
    int method_computation_cfl;

    MassMatrixParameter param_mass;

    bool same_preconditioning;
    bool first_preconditioning_computed;
    All_Preconditioner_Base<Real_wp>* ptr_real_precond;
    
    int nb_iter_per_display;
    bool display_computational_time;
    bool check_prod_stiffness_matrix;
    bool small_random_initial;
    
  public :
    Real_wp default_threshold_cfl;
    VectString param_time_scheme;
    enum { CFL_MESH_SIZE_ELT, CFL_MESH_SIZE_FACE, CFL_MESH_SIZE_VERTEX, CFL_EIGENVALUE};
    
    //! parameters of the inital condition
    VectString param_initial_condition;
    Real_wp epsilon_time;
    
    //! parameters of the time source
    Vector<string> param_time_source;

    Real_wp tinit_source; //!< initial time of the source
    Real_wp tlimit_source; //!< final time of the source
    
    //! spatial source, vector F, stored in sparse format
    Vector<Real_wp, VectSparse> sparse_vector_source, lambda_vector_source;
    Vector<Real_wp, VectSparse> dirichlet_source, dirichlet_source_unassembled;
    
    //! object containing the computation of the time source
    TimeSourceHyperbolic time_source;
    
  public :
    Real_wp GetFinalTimeSource() const;
    void SetFinalTimeSource(const Real_wp& t);
    
    VirtualTimeSource<Real_wp>* GetNewTimeSource();
    
    Vector<Real_wp, VectSparse>& GetSpaceSource();
    
    // true if f(x,t) = g(t) h(x) , false in other cases
    virtual bool SourceDoesNotDependOnTime() const;
    
    void ComputeRightHandSide();
    void FillSource(VectReal_wp& Bn) const;
            
    void SourceOnlyTime(const Real_wp& t0, int nb_deriv, Real_wp& pulse);
    void SetDirichletCondition(const Real_wp& t, int nb_deriv, VectReal_wp& b_src,
                               Real_wp alpha = 1.0);

    virtual void SetDirichletConditionSource(const Real_wp& t, int nb_deriv, VectReal_wp& b_src,
					     Real_wp alpha = 1.0) = 0;

    virtual void AddPrimitiveSourceAtTime(const Real_wp& alpha, const Real_wp& t,
					  int nb_deriv, VectReal_wp& b_src);

    virtual void AddPrimitiveTimeSource(const Real_wp& alpha, const Real_wp& t,
					int nb_deriv, VectReal_wp& b_src);

    // pure virtual functions
    virtual void AddScalarSourceAtTime(const Real_wp& alpha, const Real_wp& tn, int n, VectReal_wp& Y) = 0;
    virtual void AddVectorialSourceAtTime(const Real_wp& alpha, const Real_wp& tn, int n, VectReal_wp& Y) = 0;

    virtual void InitSource() = 0;

  private:
    VarProblem_Base& var_problem;
    VarComputationProblem& var_computation;
    DistributedProblem_Base& var_comm;
    VarOutputProblem_Base& var_output;
    VarBoundaryCondition_Base& var_boundary;
    VarSourceProblem_Base& var_source;
    
    void InitDefaultValues();

  public:
    template<class TypeEquation>
    VarInstationary_Base(HyperbolicProblem<TypeEquation>& var);

    virtual ~VarInstationary_Base();
    
    // Inline methods
    virtual bool SplitSystem() const;

    inline DistributedProblem_Base& GetDistributedProblem() { return var_comm; }
    inline const MassMatrixParameter& GetMassMatrixParameter() const { return param_mass; }

    Real_wp GetInitialTime() const;
    Real_wp GetFinalTime() const;
    
    bool IsImplicitScheme() const;
    bool FirstOrderScheme() const;
    int GetTimeSchemeType() const;
    int GetTimeSchemeOrder() const;
    
    Real_wp GetTimeStep() const;
    void SetTimeStep(const Real_wp&);
    
    MatrixVectorProductLevel& GetTimeLevelDistribution();
    const MatrixVectorProductLevel& GetTimeLevelDistribution() const;
    int GetOperatorLevel() const;

    void ForceComputationTimeStep(bool compute_dt = true);
    const string& GetFileStoringLocalDt() const;
    int GetCFL_ComputationalMethod() const;
    
    int GetNumberOfUnknowns() const;
    int GetNbScalarUnknowns() const;

    void ApplyOperatorKh(const Real_wp& alpha, const Real_wp& t, const VectReal_wp& Uh,
			 const Real_wp& beta, VectReal_wp& Yh);
    
    void ApplyOperatorDh(const Real_wp& alpha, const Real_wp& tn, const VectReal_wp& Uh,
			 const Real_wp& beta, VectReal_wp& Prod_Uh);
    
    void ApplyOperatorSh(const Real_wp& alpha, const Real_wp& tn, const VectReal_wp& Uh,
			 const Real_wp& beta, VectReal_wp& Prod_Uh);
    
    void ApplyOperatorShVectorial(const Real_wp& alpha, const Real_wp& tn, const VectReal_wp& Uh,
				  const Real_wp& beta, VectReal_wp& Prod_Uh);
    
    void ApplyOperatorDhMinusdtSh(const Real_wp& alpha, const Real_wp& tn, const VectReal_wp& Uh,
				    const Real_wp& beta, VectReal_wp& Prod_Uh);
    
    void SolveOperatorDhPlusdtSh(VectReal_wp& Y);
    void SolveOperatorDhPlusdtSh(const Real_wp&, const VectReal_wp& X, VectReal_wp& Y);
    
    void SolveOperatorDh(VectReal_wp& Y);
    void SolveOperatorDh(VectComplex_wp& Y);
    
    void SolveCholeskyDh(const SeldonTranspose& transA, VectReal_wp& Y);
    
    void SolveOperatorDh(const Real_wp&, const VectReal_wp& X, VectReal_wp& Y);
    void ApplyOperatorRhScalar(const Real_wp& alpha, const Real_wp& t, const VectReal_wp& Uh,
			       const Real_wp& beta, VectReal_wp& Vh, bool extrapol = true);
    
    void ApplyOperatorRhVectorial(const Real_wp& alpha, const Real_wp& t, const VectReal_wp& Uh,
				  const Real_wp& beta, VectReal_wp& Vh, bool extrapol = true);
    
    void ApplyOperatorBh(const Real_wp& alpha, const Real_wp& t, const VectReal_wp& Uh,
                         const Real_wp& beta, VectReal_wp& Vh);
    
    void ApplyOperatorBhMinusdtSh(const Real_wp& alpha, const Real_wp& t, const VectReal_wp& Uh,
				    const Real_wp& beta, VectReal_wp& Vh);
    
    void SolveOperatorBhPlusdtSh(VectReal_wp& Y);
    void SolveOperatorBhPlusdtSh(const Real_wp&, const VectReal_wp& X, VectReal_wp& Y);

    void SolveOperatorBh(VectReal_wp& Y);
    void SolveOperatorBh(const Real_wp&, const VectReal_wp& X, VectReal_wp& Y);

    void Assemble(VectReal_wp& Vh) const;
    
    virtual void EvaluateFunction(const Real_wp& tn, const VectReal_wp& X, VectReal_wp& Y,
				     bool invert_mass = true, bool source = true);
    
    virtual void EvaluateFunction(const Real_wp& tn, const Real_wp& alpha, int level,
				  const VectReal_wp& X, VectReal_wp& Y,
				  bool invert_mass = true, bool source = false);
    
    void GiveIterate(int nb_iter, const Real_wp& tn, VectReal_wp& Y);
    void GiveVectorialIterate(int nb_iter, const Real_wp& tn, VectReal_wp& Y);
    void GiveFinalIterate(int nb_iter, const Real_wp& tn, VectReal_wp& Y);
    void GiveNumberIterations(const Real_wp& dt, int nb_max_iter);

    // input/outputs
    void SetInputData(const string& description_field, const VectString& parameters);

    virtual void GetMemoryUsed(map<string, size_t>& var) const;
    
    virtual void WriteSnapshot(int nb_iter, const Real_wp& t, VectReal_wp& Uh,
			       bool compute_time = false);

    void LoadDataReprise(VectReal_wp& Xh);
    void SaveDataReprise(int nt, const VectReal_wp& Xh);
    
    // initialisation functions
    virtual void InitTimeIterations();
    
    int GetNbVectorialUnknowns() const;
    int GetNbVectorialUnknownsPML() const;
    
    virtual void ComputeMassMatrix();
    virtual void ComputeStiffnessMatrix();
    
    template<class MatrixSparse> void ExtractStiffnessMatrix(MatrixSparse& K);
    
    void SetLevel(int level);
    
    virtual void ComputeExtrapolationVectorial(const VectReal_wp&, const VectReal_wp&);
    virtual void ComputeExtrapolationScalar(const VectReal_wp&, const VectReal_wp&);   
    
    void ApplySecondOrderScheme();
    virtual void RunTimeIterations();
    
    // main method
    void RunAll(const string& input_file, const string& name_element, const string& name_equation,
		int num = -1);
      
    // evaluation of B U + F (we want to solve dU/dt = B U + F)
    virtual void EvaluateDerivativeFunction(const Real_wp& tn, int nb_deriv, const VectReal_wp& X,
					    VectReal_wp& Y, bool invert_mass = true, bool source = true);
    
    virtual void EvaluateDerivativeFunction(const Real_wp& tn, int nb_deriv, const Real_wp& alpha, 
					    int level, const VectReal_wp& X,
					    VectReal_wp& Y, bool invert_mass = true, bool source = false);
    
    virtual void EvaluateFunctionS(const Real_wp& tn, const VectReal_wp& X, const VectReal_wp& Xp, 
				   VectReal_wp& Y, bool invert_mass = true, bool source = true);
    
    // methods for time schemes
    void FactorizeOperatorDhPlusGammaKh(const Real_wp& alpha, const Real_wp& beta,
                                        const Real_wp& gamma);    

    void FactorizeOperatorReal(const VectReal_wp& alpha, const VectReal_wp& beta,
			       const VectReal_wp& gamma);    

    void FactorizeOperatorComplex(const VectComplex_wp& alpha, const VectComplex_wp& beta,
				  const VectComplex_wp& gamma);

    void SolveOperatorDhPlusGammaKh(const Real_wp& t, const VectReal_wp& X, VectReal_wp& Y);
    void SolveOperatorReal(const Real_wp& t, const VectReal_wp& X, VectReal_wp& Y, int n);
    
    void SolveOperatorComplex(const Real_wp& t, const VectComplex_wp& X, VectComplex_wp& Y, int);    
    
    void ApplyMassMatrix(const Real_wp& alpha, const Real_wp& t, const VectReal_wp& X,
			 const Real_wp& beta, VectReal_wp& y); 
    
    void SolveMassMatrix(VectReal_wp& X);        
    
    virtual void SetInitialVector(const Real_wp& t_begin, VectReal_wp& Y) = 0;
    
    // virtual methods
    virtual void GetMassDampingCoefficient(int i, int j, int num, int ref, Real_wp&, Real_wp&);
    virtual Real_wp GetCoefficientTauHDG(int ref) const;
    
    virtual void AddElementStiffnessTensor(int i, int& offset, const Real_wp& alpha,
					   const VectReal_wp& Uh, VectReal_wp& Vh,
                                           const Real_wp&, const Real_wp&);
    
    virtual void SolveElementStiffnessTensor(int i, int& offset, VectReal_wp& Vh, const Real_wp&);
    
    virtual void SolveAddElementStiffnessTensor(int i, int& offset, const Real_wp& alpha,
						const VectReal_wp& Uh, VectReal_wp& Vh, const Real_wp&);
    
    virtual Real_wp GetNormeSolution(const VectReal_wp& Uh) const = 0;

    virtual void WriteVectorialSnapshot(int& nb_iter, const Real_wp& t, VectReal_wp& Vh) = 0;
    virtual void CheckIdentityJacobianInPML() = 0;    
    
    virtual VirtualMassMatrix* GetNewMassMatrix() = 0;

    virtual void MltAddStiffnessScalar(const Real_wp&, int lvl, const VectReal_wp&,
				       const Real_wp&, VectReal_wp&);

    virtual void MltAddStiffnessVectorial(const Real_wp&, int lvl, const VectReal_wp&,
					  const Real_wp&, VectReal_wp&);
    
    virtual Real_wp EvaluateCFL(bool init_computation = true) = 0;
    
    virtual void ComputeCoarseFineRegion();
    virtual void ComputeLocalTimeStep(VectReal_wp&) = 0;
    virtual void DisplayCoarseFineRegion(const IVect&) = 0;

    virtual void GetExtrapolationLambdaHDG(int level, const VectReal_wp& Y, VectReal_wp& ProdL);
    virtual void MltAddStiffnessHDG(const Real_wp& alpha, int lvl, int lvl2,
				    const VectReal_wp& Y, const VectReal_wp& L,
				    const Real_wp& beta, VectReal_wp& ProdY);

    virtual void MltAddLambdaHDG(const Real_wp& alpha, int level, const VectReal_wp& L, VectReal_wp& ProdY);
    virtual void MltAddLambdaHDG(const Complex_wp& alpha, int level, const VectComplex_wp& L, VectComplex_wp& ProdY);
    
  };
  

  template<class Dimension>
  class VarInstationary_Dim : public VarInstationary_Base
  {
  private:
    DistributedProblem<Dimension>& var_problem;
    VarOutputProblem_Dim<Dimension>& var_output;
    VarSourceProblem_Cplx<Real_wp, Dimension>& var_source;
    VarBoundaryCondition_Dim<Real_wp, Dimension>& var_boundary;
    
  public:
    template<class TypeEquation>
    VarInstationary_Dim(HyperbolicProblem<TypeEquation>& var);

    inline DistributedProblem<Dimension>& GetDistributedProblem() { return var_problem; }
    
    void InitSource();
    
    // computation of temporal source
    virtual void AddScalarSourceAtTime(const Real_wp& alpha, const Real_wp& tn, int n, VectReal_wp& Y);
    virtual void AddVectorialSourceAtTime(const Real_wp& alpha, const Real_wp& tn, int n, VectReal_wp& Y);

    virtual void SetDirichletConditionSource(const Real_wp& t, int nb_deriv, VectReal_wp& b_src,
					     Real_wp alpha = 1.0);
    
    virtual void AddScalarTimeSource(const Real_wp& alpha, const Real_wp& tn, int n, VectReal_wp& Y);
    virtual void AddVectorialTimeSource(const Real_wp& alpha, const Real_wp& tn, int n, VectReal_wp& Y);

    void WriteVectorialSnapshot(int& nb_iter, const Real_wp& t, VectReal_wp& Vh);

    void CheckIdentityJacobianInPML();
    
    void SetInitialVector(const Real_wp& t_begin, VectReal_wp& Y);
    
    virtual VirtualProjectorFEM<Real_wp, Dimension>* GetNewInitialCondition() = 0;
    virtual VarInstationary_Dim<Dimension>* GetNewHyperbolicProblem() const = 0;

    virtual void ComputeLocalTimeStep(VectReal_wp&);
    virtual void DisplayCoarseFineRegion(const IVect&);
    
  };
  
  template<class Dimension>
  class VarInstationary_Fem : public VarInstationary_Dim<Dimension>
  {
  private:
    VarProblem<Dimension>& var_problem;
    
  public:
    template<class TypeEquation>
    VarInstationary_Fem(HyperbolicProblem<TypeEquation>& var);

    Real_wp GetNormeSolution(const VectReal_wp& Uh) const;
    
  };
  
  
  template<class TypeEquation>
  class VarInstationary : public VarInstationary_Fem<typename TypeEquation::TypeEquationStationary::Dimension>
  {
  protected:
    typedef typename TypeEquation::TypeEquationStationary TypeEquationStationary;    
    typedef typename TypeEquation::TypeEquationStationary::Dimension Dimension;    
    
  public:
    EllipticProblem<typename TypeEquation::TypeEquationStationary> var_harmonic;
    
  public:
    VarInstationary();
    
    HyperbolicProblem<TypeEquation>& GetLeafClass();
    
    VirtualMassMatrix* GetNewMassMatrix();
    VirtualProjectorFEM<Real_wp, Dimension>* GetNewInitialCondition();
    
    void SetInputData(const string& description_field, const VectString& parameters);

    virtual Real_wp EvaluateCFL(bool init_computation = true);

    VarInstationary_Dim<Dimension>* GetNewHyperbolicProblem() const;
    
  };
  
} // namespace Montjoie

#define MONTJOIE_FILE_VAR_INSTATIONARY_HXX
#endif
