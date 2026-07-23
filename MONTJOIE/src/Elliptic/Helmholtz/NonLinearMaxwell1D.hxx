#ifndef MONTJOIE_FILE_NON_LINEAR_MAXWELL_1D_HXX

namespace Montjoie
{
  class MaxwellProblem1D;
  
  class StiffnessOperatorMaxwell
  {
  public:
    static StiffnessOperatorMaxwell* GetNewObject(int r, MaxwellProblem1D& var);

    virtual void GetGradientMatrix() = 0;
    virtual void AddMatrixBand(int N, GlobalGenericMatrix<Real_wp>& nat_mat) = 0;
    virtual void SetFreeStorage(int N, GlobalGenericMatrix<Real_wp>& nat_mat, const Real_wp& theta) = 0;
    
    virtual void MltTiny(const Real_wp& alpha, const SeldonTranspose& trans, const VectComplex_wp& U, VectComplex_wp& V) = 0;
    virtual void MltFree(const Real_wp& alpha, const SeldonTranspose& trans, const VectComplex_wp& U, VectComplex_wp& V) = 0;

    virtual void ComputeJacobian(const VectComplex_wp& En, const VectReal_wp& scale) = 0;
    virtual void SolveTiny(VectComplex_wp& x_sol) = 0;

    virtual void MltFirstBlock(const Complex_wp& E0, const Real_wp& alpha, int offset, VectComplex_wp& b_src) = 0;

    virtual void FactorizeFast(Matrix<Real_wp, General, BandedCol>& mat_band_DF, bool hg) = 0;
    virtual void SolveFast(VectComplex_wp& Fu) = 0;

    virtual void FactorizeReal(Matrix<Real_wp, General, BandedCol>& mat_band_DF, bool hg, int i) = 0;
    virtual void SolveReal(int, VectComplex_wp& Fu) = 0;

    virtual void FactorizeComplex(Matrix<Complex_wp, General, BandedCol>& mat_band_DF, bool hg, int i) = 0;
    virtual void SolveComplex(int, VectComplex_wp& Fu) = 0;

  };

  
  template<int r>
  class TinyStiffnessOperatorMaxwell : public StiffnessOperatorMaxwell
  {
  protected:
    MaxwellProblem1D& var;
    EllipticProblem<LaplaceEquation1D>& var_laplace;
    
    //TinyBandMatrix<Real_wp, r> mat_tiny_band_stiff;
    TinyMatrix<Real_wp, General, r+1, r+1> mat_tiny_stiff;

    //TinyBandMatrix<Real_wp, r> mat_tiny_band_DF;
    TinyBlockSolver1D<Real_wp, r-1> solver_opt;

    Vector<TinyBlockSolver1D<Real_wp, r-1> > solver_real;
    Vector<TinyBlockSolver1D<Complex_wp, r-1> > solver_cplx;

    Vector<TinyBandMatrix<Complex_wp, r> > fem_harmonic_matrix, fem_harmonic_mat_lu;
    
  public:
    TinyStiffnessOperatorMaxwell(MaxwellProblem1D& var_);
    
    void GetGradientMatrix();
    void AddMatrixBand(int N, GlobalGenericMatrix<Real_wp>& nat_mat);
    void SetFreeStorage(int N, GlobalGenericMatrix<Real_wp>& nat_mat, const Real_wp& theta);

    void MltTiny(const Real_wp& alpha, const SeldonTranspose& trans, const VectComplex_wp& U, VectComplex_wp& V);
    void MltFree(const Real_wp& alpha, const SeldonTranspose& trans, const VectComplex_wp& U, VectComplex_wp& V);

    void ComputeJacobian(const VectComplex_wp& En, const VectReal_wp& scale);
    void SolveTiny(VectComplex_wp& x_sol);

    void MltFirstBlock(const Complex_wp& E0, const Real_wp& alpha, int offset, VectComplex_wp& b_src);

    void FactorizeFast(Matrix<Real_wp, General, BandedCol>& mat_band_DF, bool hg);
    void SolveFast(VectComplex_wp& Fu);

    void FactorizeReal(Matrix<Real_wp, General, BandedCol>& mat_band_DF, bool hg, int i);
    void SolveReal(int i, VectComplex_wp& Fu);

    void FactorizeComplex(Matrix<Complex_wp, General, BandedCol>& mat_band_DF, bool hg, int i);
    void SolveComplex(int i, VectComplex_wp& Fu);

  };

  
  //! main class for solving 1-D wave equation with non-linear Kerr effect
  /*!
    we consider the following equation
    eps_inf/c^2 d^2 E/dt^2 + 1/c^2 d^2 / dt^2 (\sum P_k) + \alpha dE/dt - d^2 E / dz^2
    + gamma_NL / (c^2 epsilon0) d^2/dt^2( |E|^2 E) = 0
    1/omega_k^2 d^2 P_k/dt^2 + \sigma_k dP_k/dt + P_k = \alpha_k E
    with a non-homogeneous Dirichlet condition
    E(z=0, t) = exp(-1/2 ((t-Tmax)/tau_impulse)^2 )*exp(-i omega t) 
    
    with the damping alpha = mu_0 \sigma (\sigma being the conductivity)
    and sigma_k = 1/(omega_k^2 T_k) (T_k being a damping time)
    
    This equation is discretized with high-order finite elements in z-direction,
    giving the following evolution system :
    eps_inf/c^2 Mh d^2 E/dt^2 + 1/c^2 Mh d^2 / dt^2 (\sum P_k) + Sh E + Kh E
    + gamma_NL / (c^2 epsilon0) Mh d^2/dt^2( |E|^2 E) = 0
    1/omega_k^2 d^2 P_k/dt^2 + P_k = \alpha_k E
    with Mh = mass matrix = \int \varphi_i \varphi_j dz,
    Kh = stiffness matrix = \int \nabla \varphi_i \nabla \varphi_j dz
    Sh = damping matrix = \int \sigma \varphi_i \varphi_j
    
    This evolution system can be solved with a theta-scheme :
    eps_inf/c^2 Mh (E^{n+1} - 2 E^n + E^{n-1}) / dt^2 
    + 1/c^2 Mh \sum_k (P_k^{n+1} - 2 P_k^n + P_k^{n-1}) / dt^2
    + Kh (theta E^{n+1} + (1 - 2 theta) E^n + theta E^n+1)
    + gamma_NL / (c^2 epsilon0) Mh
      (|E^{n+1}|^2 E^{n+1} - 2 |E^n|^2 E^n + |E^{n-1}|^2 E^{n-1}) / dt^2 = 0
    
    A newton method is used to solve the non-linear system in E^{n+1} at each time step
    
    The evolution system can also be solved by using a general scheme solving :
    d^2 U/dt^2 = f(t, U, dU/dt)
    The function EvaluateFunction to this aim is implemented
    The unknown U contains the displacement D and polarizations P_k

    A First-order formulation is also implemented and can be written as :
    eps_inf/c dE/dt + sigma E + 1/c d/dt (\sum P_k) - dH/dz 
    + gamma_NL / (c epsilon0) d/dt( |E|^2 E) = 0
    1/c dH/dt + sigma H - dE/dz = 0
    1/omega_k d Q_k/dt + P_k = \alpha_k E
    1/omega_k d P_k/dt - Q_k = 0

    It can be solved with a leap-frog scheme :
    
    eps_inf/c Mh (E^n+1 - E^n)/dt + Sh (E^n + E^{n+1})/2 + 1/c (\sum (P_k^(n+1) - P_k^n))/dt
    + Rh H^(n+1/2) + gamma_NL / (c epsilon0) Mh [ (|E|^2 E)^(n+1) - (|E|^2 E)^n]/dt = 0
    1/c Mh (H^(n+3/2) - H^(n+1/2))/dt + ShV (H^(n+3/2) + H^(n+1/2))/2 - Rh* E^(n+1) = 0
    
    or with any scheme (through EvaluateFunction) solving :
    dU/dt = f(t, U)
    The unknown U contains (D, P_k, H, Q_k)
    
    E is solved in complex numbers in order to take into account the different polarizations
    E = Ex + i Ey
    If the polarization is linear, E should be real
  */
  class MaxwellProblem1D : public NonLinearOpticsProblem, public VirtualOdeSystem<Complex_wp>,
			   public InputDataProblem_Base, public NonLinearEquations_Newton<Complex_wp>
  {
  public :
    //! number of degrees of freedom in z
    int nb_points_z;
    
    //! time iterates E^n, E^{n-1} and E^{n+1}
    Vector<Complex_wp> E_current, E_prev, E_next, ProdTmp;
    //! time iterates P_k^n, P_k^{n-1}, P_k^{n+1}
    Vector<Complex_wp> P_current, P_prev, P_next, rhs_polar;
    //! time iterates H^(n+1/2), H^(n+3/2)
    Vector<Complex_wp> H_current, H_next;
    Vector<Complex_wp> mass_termEn, mass_termEn_prev;
    
    //! 1-D Laplace problem
    /*!
      This object is used to construct the stiffness and mass matrices Mh, Kh
    */
    EllipticProblem<LaplaceEquation1D> var_laplace;

    //! for time-harmonic discretization
    EllipticProblem<HelmholtzEquation1D> var_helmholtz;
    
    //! stiffness matrix
    Matrix<Real_wp, General, BandedCol> mat_band_stiff;
    Matrix<Real_wp, General, ArrayRowSparse> mat_sparse_stiff;
    StiffnessOperatorMaxwell* var_stiff;
    //! type of storage for the stiffness matrix
    int type_storage_stiff;
    //! available storages
    enum{BAND_STORAGE, TINY_BAND_STORAGE, SPARSE_STORAGE, FREE_STORAGE};
    
    //! mass matrix and its inverse
    Vector<Real_wp> mat_mass, invMat_mass;
    //! damping matrix
    Vector<Real_wp> mat_sigma;

    VectReal_wp OriginalVertex; 
    
    //! jacobian matrix
    Matrix<Real_wp, General, ArrayRowSparse> mat_sparse_DF;
    Matrix<Real_wp, General, BandedCol> mat_band_DF;
    
    All_MatrixLU<Real_wp> mat_sparse_lu;
    Matrix<Real_wp> mat_first_lu;
    int type_storage_jacobian;
    //! diagonal jacobian matrix for an explicit matrix
    Vector<Real_wp> diag_inv_DF;
    //! pivot used for factorisation of jacobian matrix
    IVect pivot_DF;
    
    Vector<VectComplex_wp> mat_mass_harmonic;

    //! storage of seismogramms
    Vector<Vector<Complex_wp> > En_sismo;
    
    //! directory where output files are written
    string DOSSIER;
    
  public :
    MaxwellProblem1D();
    
    Real_wp GetTimeStep() const;
    Real_wp GetInitialTime() const;
    Real_wp GetFinalTime() const;
    int GetNbDof() const;
    
    void SetInputData(const string& description_field, const VectString& parameters);
    
    void ConstructAll(const string& input_file, int num = -1);
    
    int SetInitialCondition(bool explicit_scheme);
    inline Real_wp GetFinalTimeSource() const { return Real_wp(2)*this->Tmax; }
    
    void ConstructHarmonic();

    void SolveNewton(VectComplex_wp& source, VectComplex_wp& sol);
    void WriteDatas(VectComplex_wp& sol);
    
    Vector<Complex_wp>& GetIterate();
    
    void MltStiffness(const Real_wp& alpha, const SeldonTranspose& trans,
                      const Vector<Complex_wp>& U, Vector<Complex_wp>& V);    
    
    void EvaluateRightHandSide(const Real_wp& t);
    
    void ComputeScheme(const Vector<Complex_wp>& En, Vector<Complex_wp>& scheme);

    template<class MatrixSparse>
    void ComputeJacobianMatrix(const Vector<Complex_wp>& En, const Vector<Real_wp>& scale,
                               MatrixSparse& mat_DF);
    
    void ComputeAndFactoriseDiff(const Vector<Complex_wp>& En, const Vector<Real_wp>& scale);
    void SolveDifferential(const Vector<Complex_wp>& rhs, Vector<Complex_wp>& x_sol);
   

    // for implicit-schemes
    void FillInitialCondition(VectComplex_wp& E0);
    
    void FactorizeOperatorDhPlusGammaKh(const Real_wp& alpha, const Real_wp& beta,
                                        const Real_wp& gamma);
    void SolveOperatorDhPlusGammaKh(const Real_wp& t, const VectComplex_wp& X,VectComplex_wp & Y);

    void FactorizeOperatorReal(const VectReal_wp& alpha, const VectReal_wp& beta,
			       const VectReal_wp& gamma);

    void SolveOperatorReal(const Real_wp& t, const VectComplex_wp& X, 
			   VectComplex_wp& Y, int num_system);

    void FactorizeOperatorComplex(const VectComplex_wp& alpha, const VectComplex_wp& beta,
				  const VectComplex_wp& gamma);

    void SolveOperatorComplex(const Real_wp& t, const VectComplex_wp& X, 
			      VectComplex_wp& Y, int num_system);
    
    void SolveMassMatrix(VectComplex_wp& Y);

    void ApplyMassMatrix(const Real_wp& alpha, const Real_wp& tn, const VectComplex_wp& Uh,
			 const Real_wp& beta, VectComplex_wp& Prod_Uh);
    
    // for both implicit and explicit schemes
    void EvaluateFunction(const Real_wp& tn, const VectComplex_wp& X, VectComplex_wp& Px,
                          bool invert = true, bool source = true);
    
    void EvaluateFunctionS(const Real_wp& tn, const VectComplex_wp& X, const VectComplex_wp& Xp,
			   VectComplex_wp& ProdX, bool invert = true, bool source = true);
    

    //real source function
    void AddPrimitiveTimeSource(const Real_wp& alpha, const Real_wp& t, int nb_deriv, 
				Vector<Real_wp>& b_src);

    //complex source function
    void AddPrimitiveTimeSource(const Real_wp& alpha, const Real_wp& t, int nb_deriv, 
				Vector<Complex_wp>& b_src);

    
    void SetDirichletCondition(const Real_wp& t, int n, VectComplex_wp& Y,
			       Real_wp alpha = 1.0);
    
    Real_wp GetNorm2Vector(const Vector<Complex_wp>& x);
    
    void AdvanceScheme(int nt, const Real_wp& t, All_TimeScheme<Complex_wp>& scheme);    
    void UpdateScheme(int nt, const Real_wp& t);
    void UpdateMesh(int nt, const Real_wp& t, All_TimeScheme<Complex_wp>& scheme);
    void TranslateMesh(int N, int Ndof, int Nh, int, int, int, All_TimeScheme<Complex_wp>& scheme);
    
    template<class T>
    void TranslateIterate(Vector<T>& En, int Ndof, int Nh, int, int, int);
    
    void ForceUniquenessSismo(GridInterpolation<Dimension1>& sis, IVect&);
    
    string GetFourierOutputFile() const;
    void CloseBuffers();
    
    void ComputeRightHandSide(VectComplex_wp& source_rhs);
    
    void WriteSismo(int i);
    void WriteSnapshot(int nt, const Real_wp& t, Vector<Complex_wp>& En);

    void WriteHarmonicSismo(VectComplex_wp& Ex, VectComplex_wp& Ey, int i);

    void WriteOutputFile(const VarProblem_1D& var, Vector<Complex_wp>& En, const string& name_file);
    
  };
  
}

#define MONTJOIE_FILE_NON_LINEAR_MAXWELL_1D_HXX
#endif
