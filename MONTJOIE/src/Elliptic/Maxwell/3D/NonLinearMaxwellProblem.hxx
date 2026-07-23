#ifndef MONTJOIE_FILE_NON_LINEAR_MAXWELL_PROBLEM_HXX

namespace Montjoie
{

  //! main class for solving axisymmetric wave equation with non-linear Kerr effect
  /*!
    we consider the following equation
    eps_inf/c^2 d^2 E/dt^2 + 1/c^2 d^2 / dt^2 (\sum P_k) + \alpha dE/dt - \Delta E
    + gamma_NL / (c^2 epsilon0) d^2/dt^2( |E|^2 E) = 0
    1/omega_k^2 d^2 P_k/dt^2 + \sigma_k dP_k/dt + P_k = \alpha_k E
    with a non-homogeneous Dirichlet condition
    E(z=0, t) = exp(-1/2 ((t-Tmax)/tau_impulse)^2 )*exp(-i omega t) 
    
    with the damping alpha = mu_0 \sigma (\sigma being the conductivity)
    and sigma_k = 1/(omega_k^2 T_k) (T_k being a damping time)
    
    This equation is discretized with high-order finite element
    giving the following evolution system :
    eps_inf/c^2 Mh d^2 E/dt^2 + 1/c^2 Mh d^2 / dt^2 (\sum P_k) + Sh E + Kh E
    + gamma_NL / (c^2 epsilon0) Mh d^2/dt^2( |E|^2 E) = 0
    1/omega_k^2 d^2 P_k/dt^2 + P_k = \alpha_k E
    with Mh = mass matrix
    Kh = stiffness matrix
    Sh = damping matrix
    
    The evolution system is solved by using a general scheme solving :
    d^2 U/dt^2 = f(t, U, dU/dt)
    The function EvaluateFunction to this aim is implemented
    The unknown U contains the displacement D and polarizations P_k

    E is solved in complex numbers in order to take into account the different polarizations
    E = Ex + i Ey
    If the polarization is linear, E should be real
  */
  template<class TypeEquation>
  class NonLinearMaxwellProblem : public NonLinearOpticsProblem, public VirtualOdeSystem<Real_wp>,
				  public VirtualOdeSystem<Complex_wp>, public InputDataProblem_Base
  {
  public :
    typedef typename TypeEquation::Dimension Dimension;
    typedef typename Dimension::R_N R_N;
    typedef typename Dimension::MatrixN_N MatrixN_N;
    
    //! Laplace problem
    /*!
      This object is used to construct the stiffness and mass matrices Mh, Kh
    */
    EllipticProblem<TypeEquation> var_laplace;
    
    //! stiffness matrix
    FemMatrixFreeClass<Real_wp, TypeEquation> Glob_mat_Kh;
    //! mass matrix (diagonal) and its inverse
    Vector<Real_wp> mat_mass, invMat_mass;
    //! damping matrix (diagonal)
    Vector<Real_wp> mat_sigma;
    //! right-hand-side
    Vector<Real_wp, VectSparse> dirichlet_source;
    //! Dirichlet dof numbers
    Vector<int> dirichlet_dof;
    //! original vertices of the mesh (before translation)
    Vector<R_N> OriginalVertex; 
    //! points of the sismograms
    GridInterpolation<Dimension> grid_sismoND;
    //! adjacent processors
    IVect neighboring_proc;
    //! storage of seismogramms
    Vector<Vector<Real_wp> > En_sismo;
    int nb_components_sismo;

  public :
    NonLinearMaxwellProblem();
    
    Real_wp GetTimeStep() const;
    Real_wp GetInitialTime() const;
    Real_wp GetFinalTime() const;
    int GetNbDof() const;
    
    void SetInputData(const string& description_field, const VectString& parameters);
    
    void GenerateGridPoints(const Real_wp& xmin, const Real_wp& xmax, const Real_wp& ymin,
                       const Real_wp& ymax, int nbx, int nby, int nbz,
                            VectR2& GlobalCoord, VectReal_wp& TetaInterp);
    
    void ForceUniquenessSismo(GridInterpolation<Dimension>& grid, IVect& proc_sismo);
    
    void ComputeRightHandSide();
    
    void ConstructAll(const string& input_file, const string& name_element, int num = -1);
    
    void InitTimeIterations();
        
    template<class Vector1>
    void MltStiffness(const Vector1& U, Vector1& V);    
    
    void EvaluateFunction(const Real_wp& tn, const VectReal_wp& X, VectReal_wp& Px,
                          bool invert = true, bool source = true);
    
    void EvaluateFunctionS(const Real_wp& tn, const VectReal_wp& X, const VectReal_wp& Xp,
			   VectReal_wp& ProdX, bool invert = true, bool source = true);
    
    template<class T> Real_wp GetNorm2Vector(Vector<T>& x);
    
    template<class T>
    void AdvanceScheme(int nt, const Real_wp& t, All_TimeScheme<T>& scheme, Vector<T>&);    

    template<class T>
    void UpdateMesh(int nt, const Real_wp& t, All_TimeScheme<T>& scheme);
    
    template<class T>
    void TranslateMesh(int N, int Ndof, int Nh, int, int, int, All_TimeScheme<T>& scheme);
    
    template<class T>
    void TranslateIterate(Vector<T>& En, int Ndof, int Nh, int, int, int);
    
    void CloseBuffers();
    void WriteSismo(int i);

    template<class T> TinyVector<T, 1>
    GetInterpolate(const Vector<T>& E, int ne, const R_N& pt_loc, const R_N& pt_glob, const MatrixN_N& dfjm1);
    
    void CopyValueSismo(const TinyVector<Real_wp, 1>& val, VectReal_wp& En, int offset);
    void CopyValueSismo(const TinyVector<Complex_wp, 1>& val, VectReal_wp& En, int offset);
    void CopyValueSismo(const TinyVector<Real_wp, 3>& val, VectReal_wp& En, int offset);
    
    template<class T>
    void WriteSnapshot(int nt, const Real_wp& t, Vector<T>& En);
    
    void SetDirichletCondition(const Real_wp& t, int n, Vector<Real_wp>& Y, Real_wp alpha = 1.0);
        
  };
  
}

#define MONTJOIE_FILE_NON_LINEAR_MAXWELL_PROBLEM_HXX
#endif
