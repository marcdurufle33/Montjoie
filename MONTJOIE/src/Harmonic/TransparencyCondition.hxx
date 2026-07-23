#ifndef MONTJOIE_FILE_TRANSPARENCY_CONDITION_HXX

namespace Montjoie
{
  
  //! base class for computing radar cross section (far field)
  template<class Dimension>
  class VarComputationRCS_Base : public InputDataProblem_Base
  {
  protected :    
    //! order of integration to compute radar cross section (RCS)
    int order_integration_rcs;
    int nb_angles_RCS; //!< number of angles for RCS
    int LINE_RCS; //!< body number to compute RCS
    int type_rcs; //!< monostatic, bistatic ?
    Real_wp first_angle_RCS; //!< first angle
    Real_wp last_angle_RCS; //!< last angle
    //! if true, we compute RCS
    bool rcs_to_be_computed;
    string file_RCS; //!< file where RCS will be stored
    string file_far_field; //!< file where far field will be stored
    int ref_domain_RCS;
    int print_level;
        
    //! object used for interpolation on the surface mesh
    MeshInterpolationFEM<Dimension> var_mesh;
    
    // points outside the domain where the solution will be computed
    // through a representation formula (integral equation)
    typename Dimension::VectR_N PointsOutside;
    Real_wp deltat;
    
  public:
    enum {MONOSTATIC_RCS, BISTATIC_RCS};
    
  private:
    VarProblem<Dimension>& var_problem;

    void InitDefaultValues();

  public :
    template<class TypeEquation>
    VarComputationRCS_Base(EllipticProblem<TypeEquation>& var);

    // Inline methods
    bool RcsToBeComputed() const;
    int GetNbPointsOutside() const;
    int GetNbAngles() const;
    int GetRcsType() const;
    
    const MeshInterpolationFEM<Dimension>& GetInterpolationMesh() const;
    
    const typename Dimension::R_N& GetOutsidePoint(int i) const;
    const typename Dimension::VectR_N& GetOutsidePoint() const;
    void SetOutsidePoints(const typename Dimension::VectR_N& pts);
    void SetTimeStep(const Real_wp& dt);
    
    // to set parameters of the object
    virtual size_t GetMemorySize() const;
    void SetInputData(const string& description_field,
		      const VectString& parameters);
    
    // main method of this object
    virtual void InitComputationRCS(bool assemble = false);
    
    // almost pure virtual methods
    virtual void ComputeRCS(const VectReal_wp& U0) = 0;
    virtual void ComputeRCS(const VectComplex_wp& U0) = 0;
    
    virtual void
    ComputeIntegralRepresentation(const VectReal_wp& trace_En, const VectReal_wp& trace_Hn,
				  const MeshInterpolationFEM<Dimension>& mesh_interp,
				  const typename Dimension::R_N& pointX, VectReal_wp&) const;

    virtual void
    ComputeIntegralRepresentation(const VectComplex_wp& trace_En, const VectComplex_wp& trace_Hn,
				  const MeshInterpolationFEM<Dimension>& mesh_interp,
				  const typename Dimension::R_N& pointX, VectComplex_wp&) const;
    
  };
  
  
  //! leaf class for radar cross section (this class will be specialized)
  template<class TypeEquation>
  class VarComputationRCS : public VarComputationRCS_Base<typename TypeEquation::Dimension>
  {
  public :
    VarComputationRCS(EllipticProblem<TypeEquation>& var);
    
    virtual void ComputeRCS(const VectReal_wp& U0);
    virtual void ComputeRCS(const VectComplex_wp& U0);
    
  };  

  
  //! base class for transparent condition
  /*!
    black-box matrix for transparency condition for Helmholtz/Maxwell equation
    (A_c - A_p) X = B
    A_c : sparse matrix coming from finite element discretization
    A_p : full matrix, not stored, issued from integral representation
    matrix-vector product available for this matrix
  
    the system is solved by multiplying by A_c^{-1}
    (I - A_c^{-1} A_p) X = B
  
    the black-box matrix is A_h = (I - A_c^{-1} A_p)
    for the definition of methods, see file TransparencyCondition.cxx
  */
  class TransparencySolver_Base : public VirtualMatrix<Complex_wp>, public InputDataProblem_Base
  {
  protected :
    //! order of integration used in the integral representation
    int order_integration_transparency_condition;
        
    //! finite element solver
    All_LinearSolver& glob_solver;
    
    //! iterative solver for transparency condition
    All_IterativeSolver<Complex_wp> iterative_solver;
    
    //! if true, we use a transparency condition else a first order ABC
    bool use_transparency_condition;
    int ref_transparent_body;
    int ref_domain_transparent;
    
    //! needs to display historical of iterations
    bool display_messages;
    mutable int nb_prod;
    int first_point_to_treat, nb_points_to_treat;
    Vector<TinyVector<int, 3> > points_to_send, points_to_recv;
    GlobalGenericMatrix<Complex_wp> nat_mat;

  private :
    VarProblem_Base& var_problem;
    DistributedProblem_Base& var_comm;

    void InitDefaultValues();
    
  public :
    template<class TypeEquation>
    TransparencySolver_Base(EllipticProblem<TypeEquation>&, All_LinearSolver&);

    // Inline methods
    bool UseTransparentCondition() const;    
    void Solve(VectReal_wp & x_sol, VectReal_wp & b_source);

    void Reallocate(int, int);
    void Clear();

    void HideMessages();
    void ShowMessages();
    void ComputeSolution(VectComplex_wp& rhs, VectComplex_wp& sol) const;
    
    // other methods
    void SetInputData(const string&, const VectString&);
    
    void Solve(VectComplex_wp & x_sol, VectComplex_wp & b_source);
    
    void Init();    

    void ExchangeValues(VectComplex_wp& EnSigma, VectComplex_wp& HnSigma,
			VectComplex_wp& EnQuad, VectComplex_wp& HnQuad) const;
    
    static void ComputeGreenKernel(const R3& X, const R3& Y, const Real_wp& k, Complex_wp& phi,
			    R3_Complex_wp& grad_phi, Matrix3_3sym_Complex_wp& hessian_phi);
    
    static void ComputeGreenKernel(const R3& X, const R3& Y, const Real_wp& k, Complex_wp& phi,
				   R3_Complex_wp& grad_phi);
    
    void MltAddVector(const Real_wp& alpha, const VectReal_wp& B, const Real_wp& beta,
		      VectReal_wp& C) const;

    void MltAddVector(const Complex_wp& alpha, const VectComplex_wp& B, const Complex_wp& beta,
		      VectComplex_wp& C) const;

    void MltAddVector(const Real_wp& alpha, const SeldonTranspose&, 
                      const VectReal_wp& B, const Real_wp& beta, VectReal_wp& C) const;

    void MltAddVector(const Complex_wp& alpha, const SeldonTranspose&,
                      const VectComplex_wp& B, const Complex_wp& beta, VectComplex_wp& C) const;
    
    void MltVector(const VectReal_wp& B, VectReal_wp& C) const;
    void MltVector(const VectComplex_wp& B, VectComplex_wp& C) const;

    void MltVector(const SeldonTranspose&, const VectReal_wp& B, VectReal_wp& C) const;
    void MltVector(const SeldonTranspose&, const VectComplex_wp& B, VectComplex_wp& C) const;
    
    // pure virtual methods
    virtual size_t GetMemorySize() const = 0;
    virtual int ComputeSurfaceGammaAndAbsorbing(int ref_abc, int ref_gamma, IVect& offset_abc_proc) = 0;
    virtual void ComputeRightHandSide(const VectComplex_wp& U0, VectComplex_wp& g_source) const = 0;
    
  };

  
  template<class Dimension>
  class TransparencySolver_Dim : public TransparencySolver_Base
  {
  protected :
    typedef typename Dimension::R_N R_N;
    
    //! integration on first order absorbing boundary condition
    MeshInterpolationFEM<Dimension> interp_abc;

    //! integration for the intern boundary
    MeshInterpolationFEM<Dimension> interp_gamma;

  private :
    VarProblem<Dimension>& var_problem;
    VarSourceProblem_Base& var_source;

  public:
    template<class TypeEquation>
    TransparencySolver_Dim(EllipticProblem<TypeEquation>&, All_LinearSolver&);

    virtual size_t GetMemorySize() const;    
    int ComputeSurfaceGammaAndAbsorbing(int ref_abc, int ref_gamma, IVect& offset_abc_proc);
    
    void ComputeRightHandSide(const VectComplex_wp& U0, VectComplex_wp& g_source) const;
    
    virtual void ComputeAndStoreEnPot(const VectComplex_wp& En, const VectComplex_wp& Hn,
				      const R_N& pointX, const R_N& normaleX,
				      VectComplex_wp& EnStore, VectComplex_wp& HnStore, int k) const = 0;
    
    virtual void GetSource(const VectComplex_wp&, const VectComplex_wp&, int n,
			   const Real_wp& k, const R_N&, const R_N&, Vector<VectComplex_wp>&, int j) const = 0;
    
  };

  
  template<class Dimension, int nb_en, int nb_hn>
  class TransparencySolver_Fem : public TransparencySolver_Dim<Dimension>
  {
  protected :
    typedef typename Dimension::R_N R_N;
    
  public :
    template<class TypeEquation>
    TransparencySolver_Fem(EllipticProblem<TypeEquation>& var, All_LinearSolver& solver);
    
    virtual void ComputeAndStoreEnPot(const VectComplex_wp& En, const VectComplex_wp& Hn,
				      const R_N& pointX, const R_N& normaleX,
				      VectComplex_wp& EnStore, VectComplex_wp& HnStore, int k) const ;

    virtual void ComputeIntegralRepresentation(const VectComplex_wp&, const VectComplex_wp&,
					       const MeshInterpolationFEM<Dimension>&,
					       const R_N&, const R_N&, TinyVector<Complex_wp, nb_en>&,
					       TinyVector<Complex_wp, nb_hn>&) const = 0;
    
  };
  
  
  //! Class for transparent condition (specialized in other files)
  template<class TypeEquation>
  class TransparencySolver
    : public TransparencySolver_Fem<typename TypeEquation::Dimension,
				    TypeEquation::nb_components_en, TypeEquation::nb_components_hn>
  {    
    typedef typename TypeEquation::Dimension Dimension;
    typedef typename Dimension::R_N R_N;
    enum{nb_en = TypeEquation::nb_components_en, nb_hn = TypeEquation::nb_components_hn};
    
  public :
    
    TransparencySolver(EllipticProblem<TypeEquation>& var, All_LinearSolver& solver);
    
    void GetSource(const VectComplex_wp&, const VectComplex_wp&, int n,
		   const Real_wp& k, const R_N&, const R_N&, Vector<VectComplex_wp>&, int j) const;

    void ComputeIntegralRepresentation(const VectComplex_wp&, const VectComplex_wp&,
				       const MeshInterpolationFEM<Dimension>&,
				       const R_N&, const R_N&, TinyVector<Complex_wp, nb_en>&,
				       TinyVector<Complex_wp, nb_hn>&) const;
    
  };  
  
} // namespace Montjoie

#define MONTJOIE_FILE_TRANSPARENCY_CONDITION_HXX
#endif

