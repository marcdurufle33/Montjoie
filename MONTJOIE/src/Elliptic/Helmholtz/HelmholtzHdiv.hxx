#ifndef MONTJOIE_FILE_HELMHOLTZ_HDIV_HXX

namespace Montjoie
{

  //! base class for resolution of Helmholtz equation with H(div) formulation
  template<class T, class Dim>
  class HelmholtzEquationHdiv_Base : public GenericEquation<T>
  {
  public :
    typedef Dim Dimension;
    
    static bool store_dfjm1;
    
    enum {nb_unknowns = 1, nb_unknowns_scal = 1, nb_unknowns_hdg=1,
	  nb_components_en = 1, nb_components_hn = 1, nb_unknowns_vec = 1,
          type_element = 3};
    
    static inline bool SymmetricGlobalMatrix() { return true; }
    static inline bool SymmetricElementaryMatrix() { return false; }
    
    static inline bool ComputeDFjm1() { return store_dfjm1; }
    
    template<class TypeEquation>
    static void ComputeMassMatrix(EllipticProblem<TypeEquation>& var,
                                  int i, const ElementReference_Dim<Dimension>&);
    
  };

  template<class T, class Dimension>
  bool HelmholtzEquationHdiv_Base<T, Dimension>::store_dfjm1(false);


  //! resolution of Helmholtz equation with H(div) formulation
  template<class Dimension>
  class HelmholtzEquationHdiv : public HelmholtzEquationHdiv_Base<Complex_wp, Dimension>
  {
  public :
  };


  //! resolution of Laplace equation with H(div) formulation
  template<class Dimension>
  class LaplaceEquationHdiv : public HelmholtzEquationHdiv_Base<Real_wp, Dimension>
  {
  public :
  };

  
  //! Base class for indexes in H(div) formulation of Helmholtz equation
  template<class Dimension>
  class VarHelmholtzHdiv_Base
  {
    typedef typename Dimension::R_N R_N; //!< point in R2 or R3
    typedef typename Dimension::VectR_N VectR_N; //!< vector of points
    typedef typename Dimension::MatrixN_Nsym MatrixN_Nsym; //!< 2x2 or 3x3 matrix

  private:
    VarProblem<Dimension>& var_problem;
    
  public:
    //! damping
    Vector<ScalarPhysicalIndice<Dimension, Real_wp> > ref_sigma;
    //! density rho
    Vector<ScalarPhysicalIndice<Dimension, Real_wp> > ref_rho0;
    //! sound speed and pressure
    Vector<ScalarPhysicalIndice<Dimension, Real_wp> > ref_c0, ref_p0;
    Vector<VectorPhysicalIndice<Dimension, Dimension::dim_N, Real_wp> > ref_grad_p0;
    Vector<VectReal_wp> eval_rho, eval_rhoC2, eval_sigma;
    Vector<VectR_N> grad_p0;
    Vector<Vector<MatrixN_Nsym> > hessian_p0;

  public:
    template<class TypeEquation>
    VarHelmholtzHdiv_Base(EllipticProblem<TypeEquation>& var);    

    // physical indices to initialize
    void InitIndices(int n);
    void SetIndices(int i, const VectString& parameters);
    void SetPhysicalIndex(const string& name_media, int i, const VectString& parameters);
    string GetPhysicalIndexName(int m) const;
    int GetNbPhysicalIndices() const;
    bool IsVaryingMedia(int i) const;
    Real_wp GetVelocityOfMedia(int ref) const;
    Real_wp GetVelocityOfInfinity() const;
    
    void GetVaryingIndices(Vector<PhysicalVaryingMedia<Dimension, Real_wp>* >&  rho_real, IVect& num_ref,
			   IVect& num_index, IVect& num_component, Vector<bool>& compute_grad,
			   Vector<bool>& compute_hess);
    
    void ComputePhysicalCoefficients();

  };
  
  
  //! Base class for the solution of Helmholtz equation with H(div) formulation
  /*!
    The considered equation is the following one :
    -\rho \omega^2 u - div( \rho c^2 \nabla u) + (div u) \nabla p + Hess(p) u - \nabla( u \cdot \nabla p) = f
    rho, c, p are background density sound velocity and pressure
    u is the unknown displacement
    
    Adimensionalization :
    omega = omega_0 Omega, x = x_0 X, p = p_0 P, rho = rho_0 Rho, c = c_0 C
    du/dx = 1/x_0 du/dX
    
    - (\rho_0 \omega_0^2 x_0^2 / p_0) Rho Omega^2 u - (\rho_0 c_0^2)/p_0 div(Rho C^2 \nabla) 
    + (div u) \nabla P + Hess(P) u - \nabla( u \cdot \nabla P) = f/p_0

    As a result, we have rho_0 c_0^2 / p_0 = 1 => rho_0 / p_0 = 1/c_0^2
    and omega_0 = c_0 x_0
   */
  template<class TypeEquation>
  class VarHelmholtzHdiv_Eq : public VarHelmholtzHdiv_Base<typename TypeEquation::Dimension>,
			      public VarHarmonic<TypeEquation>
  {
  public :
    typedef typename TypeEquation::Dimension Dimension;
    typedef typename Dimension::R_N R_N; //!< point in R2 or R3
    typedef typename Dimension::VectR_N VectR_N; //!< vector of points
    typedef typename Dimension::MatrixN_Nsym MatrixN_Nsym; //!< 2x2 or 3x3 matrix
    typedef typename TypeEquation::Complexe Complexe; //!< real or complex
    
    
    Vector<Vector<Real_wp> > Glob_matMass_Dh;
    Vector<Vector<MatrixN_Nsym> > Glob_matMass_Bh;
    Vector<Vector<TinyMatrix<Complexe, Symmetric, Dimension::dim_N, Dimension::dim_N> > >
    Glob_matMass_AhSigma, Glob_matMass_Ah;
    Vector<Vector<R_N> > Glob_matMass_Vh;
    
    VarHelmholtzHdiv_Eq();

    void InitIndices(int n);
    void SetIndices(int i, const VectString& parameters);
    void SetPhysicalIndex(const string& name_media, int i, const VectString& parameters);
    string GetPhysicalIndexName(int m) const;
    int GetNbPhysicalIndices() const;
    bool IsVaryingMedia(int i) const;
    bool IsVaryingMedia(int m, int i) const;
    Real_wp GetVelocityOfMedia(int ref) const;
    Real_wp GetVelocityOfInfinity() const;
    
    void GetVaryingIndices(Vector<PhysicalVaryingMedia<Dimension, Complex_wp>* >&  rho_cplx,
			   Vector<PhysicalVaryingMedia<Dimension, Real_wp>* >&  rho_real, IVect& num_ref,
			   IVect& num_index, IVect& num_component, Vector<bool>& compute_grad,
			   Vector<bool>& compute_hess);
    
    void ComputePhysicalCoefficients();
    int GetBoundaryConditionId(const IVect& ref, int pos, const VectString&, bool& periodic);

    template<class T>
    void ComputeElementaryMatrixHelm(int iquad, IVect& num_dof, VirtualMatrix<T>& mat_interac,
				     const GlobalGenericMatrix<T>& nat_mat,
				     const ElementReference<Dimension, 3>& Fb);
    
    void AllocateMassMatrices();

  };


  //! source for the scattering of an incident wave
  /*!
    the diffracted field is computed
  */
  template<class Dimension>
  class DiffractedWaveSource_HelmHdiv
    : public VirtualSourceFEM<Complex_wp, Dimension>
  {
  protected :
    typedef typename Dimension::R_N R_N; //!< point in R2 or R3
    typedef typename Dimension::VectR_N VectR_N; //!< vector of points
    typedef typename Dimension::MatrixN_Nsym MatrixN_Nsym; //!< 2x2 or 3x3 symmetric matrix
    typedef typename Dimension::R_N_Complex_wp R_N_Complex_wp; //!< "complex" point
    typedef typename Dimension::MatrixN_Nsym_Complex_wp MatrixN_Nsym_Complex_wp;
    //!< 2x2 or 3x3 complex symmetric matrix
    
    //! coef_vol -> coefficient for volumic integrals  \f$ \int f \varphi \f$
    Complex_wp coef_vol;
    // physical indices : rho, mu, rho0, mu0, sigma_rho, sigma_mu
    Complex_wp invRho0; //!< rho in empty media
    Complex_wp invRho; //!< rho in element
    Complex_wp invMu0; //!< mu in empty media
    Complex_wp omega2; //!< omega^2
    MatrixN_Nsym_Complex_wp invMu; //!< mu in element
    
    //! coefficient for volumic integrals coef_grad \f$\int f \cdot grad(phi) \f$
    MatrixN_Nsym_Complex_wp coef_grad;
    
    // in order to treat high conductivity boundary conditions
    int reference_condition; //!< boundary condition
    Complex_wp coef_invMu0; //!< coefficient
    Complex_wp coef_invMu; //!< coefficient
    bool varying_rho; //!< variable rho ?
    bool varying_mu; //!< variable mu ?
    
  private:
    IncidentWaveField<Complex_wp, Dimension>& incident_wave;
    const VarProblem<Dimension>& var_problem;
    
    void InitDefaultValues();
    
  public :
    template<class TypeEquation>
    DiffractedWaveSource_HelmHdiv(const EllipticProblem<TypeEquation>& vars_helm,
				  IncidentWaveField<Complex_wp, Dimension>&);
    
    void InitElement(int num_elem, const VectR_N& s);
    
    // dirichlet condition
    void EvaluateFunction(int i, int j, const R_N& x, VectComplex_wp& f);
    
    // volumic integral \int f \phi
    void EvaluateVolumetricSource(int i, int j, const R_N& x, VectComplex_wp& f);
    bool IsNonNullVolumetricSource(const VectR_N& s);
    
    // volumic integral \int f grad(\phi)
    void EvaluateGradientSource(int i, int j, const R_N& x, VectComplex_wp& f);
    bool IsNonNullGradientSource(const VectR_N& s);
    
    // initialization for surfacic integration
    void InitSurface(int i, int num_face, int num_elem, int num_loc);
    
    // surfacic integral \int_\Gamma f \phi
    void EvaluateSurfacicSource(int k, const SetPoints<Dimension>& PointsElem,
				 const SetMatrices<Dimension>& MatricesElem, VectComplex_wp& f);
    bool IsNonNullSurfacicSource(int cond_ref);
    
    // surfacic integral \int_\Gamma f \cdot \nabla_\Gamma \phi
    void EvaluateSurfacicSourceGradient(int k,const SetPoints<Dimension>& PointsElem,
					 const SetMatrices<Dimension>& MatricesElem, VectComplex_wp& f);
    
    bool IsNonNullSurfacicSourceGradient(int cond_ref);
    
  };


  //! projector for incident wave and DG
  template<class Dimension>
  class IncidentWaveProjector_HelmHdiv : public IncidentWaveProjector<Complex_wp, Dimension>
  {
  private:
    typedef typename Dimension::R_N R_N;

    const VarHelmholtzHdiv_Base<Dimension>& var_helm;
    const VarProblem<Dimension>& var_problem;

  public :
    template<class TypeEquation>
    IncidentWaveProjector_HelmHdiv(const EllipticProblem<TypeEquation>& var,
				   IncidentWaveField<Complex_wp, Dimension>& u_inc);
    
    void EvaluateFunction(int i, int j, const R_N& x, VectComplex_wp& f);
    
  };


  //! source for the scattering of an incident wave
  template<class Dimension>
  class DiffractedWaveSource<HelmholtzEquationHdiv<Dimension> >
    : public DiffractedWaveSource_HelmHdiv<Dimension>
  {
  public :
    template<class TypeEquation>
    DiffractedWaveSource(const EllipticProblem<TypeEquation>& vars_helm,
			 IncidentWaveField<Complex_wp, Dimension>& u_inc);
    
  };
  

  //! base class for impedance boundary condition
  template<class Dimension>
  class ImpedanceABC<Complex_wp, HelmholtzEquationHdiv<Dimension> >
    : public ImpedanceFunction_Base<Complex_wp, Dimension>
  {
  private:
    const VarHelmholtzHdiv_Base<Dimension>& var_helm;
    const VarProblem<Dimension>& var_problem;
    int order; 
    Real_wp delta;
  public :
    template<class TypeEquation>
    ImpedanceABC(const EllipticProblem<TypeEquation>& var);
    
    void EvaluateImpedancePhi(int i, int num_elem, int num_edge, int num_loc, int k,
			      const GlobalGenericMatrix<Complex_wp>& nat_mat, int ref,
			      const SetPoints<Dimension>& Pts, const SetMatrices<Dimension>& Mat);
    
  };

  //! class for the solution of Helmholtz equation with H(div) elements
  template<class Dimension>
  class EllipticProblem<HelmholtzEquationHdiv<Dimension> >
    : public VarHelmholtzHdiv_Eq<HelmholtzEquationHdiv<Dimension> >
  {
  public:
    IncidentWaveProjector<Complex_wp, Dimension>*
    GetNewIncidentProjector(int n, const Vector<VectString>& param,
			    IncidentWaveField<Complex_wp, Dimension>& u_inc) const;

    void ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Real_wp>& mat_elem,
				 CondensationBlockSolver_Base<Real_wp>&,
				 const GlobalGenericMatrix<Real_wp>& nat_mat);

    void ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Complex_wp>& mat_elem,
				 CondensationBlockSolver_Base<Complex_wp>&,
				 const GlobalGenericMatrix<Complex_wp>& nat_mat);
    
  };

  //! class for the solution of Laplace equation with H(div) elements
  template<class Dimension>
  class EllipticProblem<LaplaceEquationHdiv<Dimension> >
    : public VarHelmholtzHdiv_Eq<LaplaceEquationHdiv<Dimension> >
  {
  public:
    void ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Real_wp>& mat_elem,
				 CondensationBlockSolver_Base<Real_wp>&,
				 const GlobalGenericMatrix<Real_wp>& nat_mat);

    void ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Complex_wp>& mat_elem,
				 CondensationBlockSolver_Base<Complex_wp>&,
				 const GlobalGenericMatrix<Complex_wp>& nat_mat);
    
  };

}

#define MONTJOIE_FILE_HELMHOLTZ_HDIV_HXX
#endif
