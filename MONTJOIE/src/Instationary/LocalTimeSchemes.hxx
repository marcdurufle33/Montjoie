#ifndef MONTJOIE_FILE_LOCAL_TIME_SCHEMES_HXX

namespace Montjoie
{

  //! base class for local time schemes
  class LocalTimeScheme_Base
  {
  protected:
    bool local_implicit;
    int type_splitting;
    Real_wp dt_splitting;
    string file_with_number_element;
    int p; // number of local time step
    IVect num_implicit;
    
  public :
    enum {FILE_SPLITTING, AUTO_SPLITTING};

    LocalTimeScheme_Base();

    // inline methods
    int GetSplittingType() const;
    Real_wp GetSplittingDt() const;

    string GetNumberElementFileName() const;
    void SetImplicitElement(const IVect& num);

    // other methods
    void SetInputSplitting(const Vector<string>& param, int& nb);
    
  };


  //! base class for local ode scheme
  class LocalOdeScheme_Base
  {
  protected:
    Vector<int> dof_local_L, dof_local_uv;
    Vector<int> dof_proche_uv;
    int nb_dof_uv_local, nb_dof_L_local;
    int nb_dof_uv_proche;
    Vector<int> num_elt_local;

    Vector<int> num_elt_extended;
    int nb_dof_uv_extended, nb_dof_L_extended;
    Vector<int> dof_extended_L, dof_extended_uv;

    int nb_dof_uv_extended_nonloc;
    Vector<int> dof_extended_uv_nonloc;
    
    bool matrix_stored;
    int num_level_vol, num_level_ch;

    DistributedMatrix<Real_wp, General, ArrayRowSparse> A12_sub, A21_sub, A22_sub;
    
    VectReal_wp& x_glob, &y_glob, &L_glob, &F_glob;

    Vector<int> ProcNumberLocalL;   
    Vector<IVect> RowNumberLocalL;
    DistributedMatrixIntegerArray mat_direct_info; 
    
  public:
    LocalOdeScheme_Base(VectReal_wp&x , VectReal_wp& y, VectReal_wp& L, VectReal_wp& F);

    // inline methods
    int GetNbLocalDof() const;
    int GetNbLocalDofLambda() const;
    int GetNbDofProche() const;
    int GetVolumeLevel() const;
    
    const Vector<int>& GetLocalDof() const;
    const Vector<int>& GetDofProche() const;
    
    // other methods;
    void CreateListIndex(int& nb_dof, IVect& liste_dof, IVect& IndexDof);
    int FindCommonDofs(const Vector<int>& dof_far, const Vector<int>& dof_loc,
		       Vector<int>& dof_common_fine, Vector<int>& dof_common_coarse);
    
    void Init(const Vector<int>& num_element, VirtualOdeSystem<Real_wp>& sys,
	      Vector<IVect>& lvl_vol, Vector<IVect>& lvl_surface,
	      bool compute_direct_info = false);

    void AssembleLocal(VirtualOdeSystem<Real_wp>& sys, VectReal_wp& L);
    void AssembleLocal(VirtualOdeSystem<Real_wp>& sys, VectComplex_wp& L);
    
    void ApplyOperatorAP(const VectReal_wp& X, VectReal_wp& Y,
			 VirtualOdeSystem<Real_wp>& sys, bool invert_mass);

    template<class T>
    void SolveOperatorM(VirtualOdeSystem<Real_wp>& sys, Vector<T>& Y,
			const Vector<int>& liste_dof, bool global_vec);

    void AddOperatorM(VirtualOdeSystem<Real_wp>& sys, const Real_wp& alpha,
		      const VectReal_wp& X, VectReal_wp& Y,
		      const Vector<int>& liste_dof);
    
  };


  //! local problem
  class VarComputationProblemLocal : public VarComputationProblem_Base
  {
    VarComputationProblem& var_global;
    LocalOdeScheme_Base& var_local;
    MatrixVectorProductLevel& level_vol;
    Vector<int>& IndexProche, &IndexLambda;
    MatrixVectorProductLevel level_inside;
    
  public :
    VarComputationProblemLocal(VarComputationProblem& var_global_,
			       LocalOdeScheme_Base& var_local_,
			       MatrixVectorProductLevel& level_vol_,
			       IVect& IndexProche, IVect& IndexLambda);

    int GetNbElt() const;
    int GetNbRows() const;
    int GetPrintLevel() const;
    
    void ComputeElementaryMatrix(int, IVect&, VirtualMatrix<Real_wp>&, CondensationBlockSolver_Base<Real_wp>&,
				 const GlobalGenericMatrix<Real_wp>&);
    
    void ComputeElementaryMatrix(int, IVect&, VirtualMatrix<Complex_wp>&, CondensationBlockSolver_Base<Complex_wp>&,
				 const GlobalGenericMatrix<Complex_wp>&);
    
    void GetInternalNodesElement(int i, int nb_dof_loc, int& nb_dof_edges, int& nb_dof_int,
				 Vector<int>& intern_node) const;
    
  };
    

  //! class for specifying a local ode scheme
  class LocalOdeSchemeHDG : public LocalOdeScheme_Base, public VirtualOdeSystem<Real_wp>
  {
    Globatto<Real_wp> lob;
    Real_wp dt, tn;
    bool presence_source;
    VirtualOdeSystem<Real_wp>& global_ode;
    const Vector<VectReal_wp>* wp_i_ptr;
    const Vector<Vector<Real_wp, VectSparse> >* Fp_i_ptr;
    const Vector<int>* dof_wp_ptr;
    const Vector<int>& dof_local;
    VectReal_wp F, y_imp;

    bool static_condensation;
    
    CondensationBlockSolver_Base<Real_wp>* cond_real_solver;    
    GlobalGenericMatrix<Real_wp> nat_mat_real;
    All_MatrixLU<Real_wp> mat_lu_real;

    Vector<CondensationBlockSolver_Base<Complex_wp>* > cond_cplx_solver;
    Vector<All_MatrixLU<Complex_wp> > mat_lu_cplx;
    Vector<GlobalGenericMatrix<Complex_wp> > nat_mat_cplx;

    Vector<Complex_wp> F_cplx, x_cplx, L_cplx;
    
  public:
    LocalOdeSchemeHDG(const VectReal_wp& ci, const Real_wp& dt0, VirtualOdeSystem<Real_wp>&,
		      const Vector<int>& dof_l, VectReal_wp& x , VectReal_wp& y, VectReal_wp& L,
		      VectReal_wp& Fg);

    virtual ~LocalOdeSchemeHDG();
    
    void EvaluateFunction(const Real_wp& t, const VectReal_wp& y, VectReal_wp& f,
			  bool invert_mass = true, bool source = true);

    void SetDirichletCondition(const Real_wp& t, int n, VectReal_wp& Y, Real_wp alpha = 1.0);
    Real_wp GetFinalTimeSource() const;

    void SetSource(const Real_wp& t, const Vector<VectReal_wp>& wp_i, const Vector<int>& dof_wp,
		   bool presence_source_, const Vector<Vector<Real_wp, VectSparse> >& Fp_i);

    void SolveMassMatrix(Vector<Real_wp>& X);
    void ApplyMassMatrix(const Real_wp& alpha, const Real_wp& t,
			 const Vector<Real_wp>& x, const Real_wp& beta, Vector<Real_wp>& y);
    
    void AddPrimitiveTimeSource(const Real_wp& alpha, const Real_wp& t,
				int nb_deriv, Vector<Real_wp>& b_src);

    template<class T>
    void FactorizeGeneric(const T& a, const T& b, const T& c,
			  CondensationBlockSolver_Base<T>*& cond_solver,
			  GlobalGenericMatrix<T>& nat_mat, All_MatrixLU<T>& mat_lu);
    
    void FactorizeOperatorDhPlusGammaKh(const Real_wp& a, const Real_wp& b, const Real_wp& c);

    template<class T>
    void SolveGeneric(CondensationBlockSolver_Base<T>*& cond_solver,
		      GlobalGenericMatrix<T>& nat_mat, All_MatrixLU<T>& mat_lu, Vector<T>& x,
		      Vector<T>& F_g, Vector<T>& L_g, Vector<T>& x_g);
    
    void SolveOperatorDhPlusGammaKh(const Real_wp& t, const VectReal_wp& f, VectReal_wp& y);

    void FactorizeOperatorComplex(const VectComplex_wp& alpha, const VectComplex_wp& beta,
				  const VectComplex_wp& gamma);

    void SolveOperatorComplex(const Real_wp& t, const VectComplex_wp& X, VectComplex_wp& Y,
			      int num_system);

    void GetMemoryUsed(map<string, size_t>& var) const;
    
  };
  

  //! local time-stepping scheme as proposed by Grote
  template<class T>
  class LocalTimeSteppingRK : public LocalTimeScheme_Base,
			      public VirtualTimeScheme<T>
  {
  protected:
    Real_wp dt;
    LinearRungeKutta_Iterator<T> explicit_scheme;
    VirtualTimeScheme<T>* implicit_scheme;
    Vector<T> Un;
    LocalOdeSchemeHDG* fine_ode;
    LocalOdeScheme_Base coarse_ode;
    Vector<IVect> num_level_volume, num_level_surface;
    IVect num_explicit;

    Matrix<Real_wp> DerMat;
    
    //Vector<VectReal_wp> EvalF, EvalFproche;
    Vector<Vector<Real_wp, VectSparse> > EvalF, EvalFproche;
    Vector<VectReal_wp> w_n_i, wp_i;
    VectReal_wp w_tild, y_imp, y_far, AP_y, AImP_y;

    VectReal_wp x_glob, y_glob, L_glob, F_glob;
    Vector<int> dof_common_fine, dof_common_coarse, dof_original_coarse, dof_fine_among_extended;
    
  public:
    LocalTimeSteppingRK();
    ~LocalTimeSteppingRK();

    UnivariatePolynomial<Real_wp> GetStabilityFunction() const;
    
    void SetInputData(int r, const Vector<string>& param);
    
    void SetInitialCondition(const Real_wp& t, const Real_wp& dt0,
			     Vector<T>& U0, VirtualOdeSystem<T>& sys);
    
    void Advance(const Real_wp& t, int n, VirtualOdeSystem<T>& sys);
    
    void ChangeTimeStep(const Real_wp&);
    
    Vector<T>& GetIterate();
    const Vector<T>& GetIterate() const;

    int GetNumberOfIterates() const;
    Vector<T>& GetIterate(int k);

    void Clear();
    void ClearFirst(Vector<T>& Un);

    size_t GetMemorySize() const;
    void GetMemoryUsed(map<string, size_t>& var) const;
    
  };


  template<class T>
  class LocalImexRK : public LocalTimeScheme_Base,
                      public VirtualTimeScheme<T>
  {
  protected:
    Real_wp dt;

    LocalOdeSchemeHDG* fine_ode;
    LocalOdeScheme_Base coarse_ode;

    VectReal_wp bi, ci;
    Matrix<Real_wp> Aexp, Aimp;
    
    Vector<T> Yn, ui, Fi, u_fine, y_imp, y_far;
    Vector<Vector<T> > AImP_y, AP_y;
    
    VectReal_wp x_glob, y_glob, L_glob, F_glob;
    Vector<int> dof_common_fine, dof_common_coarse, dof_original_coarse, dof_fine_among_extended;
    
  public:
    LocalImexRK();
    ~LocalImexRK();
    
    Vector<IVect> num_level_volume, num_level_surface;
    IVect num_explicit;

    UnivariatePolynomial<Real_wp> GetStabilityFunction() const;
    
    void SetInputData(int r, const Vector<string>& param);
    
    void SetInitialCondition(const Real_wp& t, const Real_wp& dt0,
                             Vector<T>& U0, VirtualOdeSystem<T>& sys);
    
    void Advance(const Real_wp& t, int n, VirtualOdeSystem<T>& sys);

    void ChangeTimeStep(const Real_wp& dt0);

    Vector<T>& GetIterate();
    const Vector<T>& GetIterate() const;

    int GetNumberOfIterates() const;
    Vector<T>& GetIterate(int k);

    void Clear();
    void ClearFirst(Vector<T>& Un);

    size_t GetMemorySize() const;
    void GetMemoryUsed(map<string, size_t>& var) const;
    
  };

  /*
  //! local time-stepping scheme proposed by Imperiale
  template<class T>
  class LocalImperialeScheme_Iterator
  {
  public :
    VectReal_wp Yn, Ynm1, Fy, Py, ImPy, W, Q0, Q1;
    Vector<VectReal_wp> Qn;
    Real_wp dt, safety_coef;
    IVect level;
    Vector<VectReal_wp> coefT;
    Vector<VectReal_wp> coefP;
    bool two_level, change_dt;
    Real_wp new_dt, dt_opt;
    bool unstable_algo;
    int order;
    Real_wp alpha1, alpha2;
    bool rk_scheme;
    DistributedMatrix<Real_wp, General, ArrayRowSparse> Glob_mat_Sh;

    
    LocalImperialeScheme_Iterator();
    
    void SetInputData(const VectString& param);
    
    void SetInitialCondition(const T& t, const T& deltat, const Vector<T>& local_dt,
                             Vector<T>& Y0, VirtualOdeSystem<T>& sys);
    
    // main method
    void Advance(const T& t, int n, VirtualOdeSystem<T>& sys);
    
    void ApplyOperatorP(int k, VectReal_wp& Y, VirtualOdeSystem<T>& sys);
    
  };
  

  template<class T>
  void RunTimeScheme(const T& t0, const T& tf, T& dt, Vector<T>& local_deltat,
		     Vector<T>& Y0, VirtualOdeSystem<T>& sys,
		     LocalImperialeScheme_Iterator<T>& local_time_scheme);


  //! local time-stepping scheme proposed by Piperno
  template<class T>
  class LocalPipernoScheme_Iterator
  {
  public :
    VectReal_wp En, Hn, ProdEn, ProdHn;
    Real_wp new_dt, safety_coef;
    IVect level;
    bool change_dt;
    int level_max;
    IVect NbEltPerLevel;
    
    LocalPipernoScheme_Iterator();
    
    void SetInputData(const VectString& param);
    
    void SetInitialCondition(const T& t, const T& deltat, const Vector<T>& local_dt,
                             Vector<T>& Y0, Vector<T>&, VirtualOdeSystem<T>& sys);
    
    void Advance(const T& t, int n, VirtualOdeSystem<T>& sys);

    void AdvanceLevel(const T& t, const Real_wp& dt, int level, VirtualOdeSystem<T>& sys);
    
  };

  template<class T>
  void RunTimeScheme(const T& t0, const T& tf, T& dt, Vector<T>& local_deltat,
		     Vector<T>& E0, Vector<T>& H0, VirtualOdeSystem<T>& sys,
		     LocalPipernoScheme_Iterator<T>& local_time_scheme);
  */
}

#define MONTJOIE_FILE_LOCAL_TIME_SCHEMES_HXX
#endif
