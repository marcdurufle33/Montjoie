#ifndef MONTJOIE_FILE_HELMHOLTZ_1D_HXX

namespace Montjoie
{
  
  //! class to specify the resolution of Helmholtz equation in 1-D
  class HelmholtzEquation1D : public GenericEquation<Complex_wp>
  {
  public :
    enum {nb_unknowns = 1};
    
    static bool ComputeDFjm1();
    
  };


  //! class to specify the resolution of Laplace equation in 1-D
  class LaplaceEquation1D : public GenericEquation<Real_wp>
  {
  public :
    enum {nb_unknowns = 1};
    
    static bool ComputeDFjm1();
    
  };

    
  //! class for solving Helmholtz (or Laplace) equation in 1-D  
  /*!
    For second-order formulation, we consider the following equation :
    (coef_rho + coef_sigma) u - d/dx ( mu du/dx ) = f
    with coef = -omega^2 rho for Helmholtz equation, and coef = rho for Laplace equation
    coef_sigma = -i omega sigma for Helmholtz equation and coef_sigma = sigma for Laplace equation
    
    The associated evolution problem is equal to :
    coef_rho d^2 u/dt^2 + coef_sigma du/dt  - d/dx( mu du/dx ) = f
    
    For first-order formulation, we consider the following set of equations
    (coef_u + coef_sigma_u) u - d/dx(v) = f
    (coef_v + coef_sigma_v) v - d/dx(u) = 0
    with coef_u = -i omega rho, coef_v = -i omega mu^{-1} for Helmholtz equation
    and coef_u = rho, coef_v = mu^{-1} for Laplace equation
    coef_sigma_u = sigma and coef_sigma_v = 0 except inside PML elements
  */
  template<class Complexe>
  class VarHelmholtz_1D : public VarProblem_1D
  {
  public :
    typedef Vector<Complexe> VectComplexe; //!< "complex" vector
    
    //! physical index in empty media
    Complexe rho0, mu0, sigma0;
    Complexe kwave;
    //! physical properties
    Vector<PhysicalVaryingMedia<Dimension1, Complexe> > ref_rho, ref_mu, ref_sigma, ref_true_rho;
    //! impedance on the extremities
    Complexe beta_impedance_left, beta_impedance_right;
    bool helio_scaling_model;
    Vector<bool> log_rho_read;
    //! coefficients J_i mu \omega_k
    Vector<VectComplexe> mat_mass_mu; 
    Vector<VectComplexe> mat_mass, mat_flow, mat_massV, mat_mass_radius;
    //!< coefficients \omega_k coef_u J_i and \omega_k coef_v J_i
    Vector<VectComplexe> mat_sigma, mat_sigmaV;
    //!< coefficients \omega_k coef_sigma_u J_i and \omega_k coef_sigma_v Ji

    Complexe left_dirichlet_condition, right_dirichlet_condition;
       
  public :
    VarHelmholtz_1D();
    
    // Inline methods
    void GetIkwave(Real_wp& ikwave);
    void GetIkwave(Complex_wp& ikwave);

    void ComputeElementaryMatrix(int i, IVect& num, VirtualMatrix<Real_wp>& A,
				 CondensationBlockSolver_Base<Real_wp>&,
				 const GlobalGenericMatrix<Real_wp>& nat_mat);

    void ComputeElementaryMatrix(int i, IVect& num, VirtualMatrix<Complex_wp>& A,
				 CondensationBlockSolver_Base<Complex_wp>&,
				 const GlobalGenericMatrix<Complex_wp>& nat_mat);

    // Other methods
    void SetInputData(const string& description_field, const VectString& parameters);
    void GetPhysicalCoefficient(Complexe& rho, Complexe& mu, Complexe& sigma,
                                Complexe& v0, int i, int j, int ref);
    
    void ComputePhysicalCoefficients();
    void ComputeMassMatrix();
    
    void GetMassMatrix(Vector<Complexe>& Dh);
    void GetDampingMatrix(Vector<Complexe>& Dh);

  protected :
    template<class T>
    void ComputeElementaryMatrixGen(int iquad, IVect& num_dof,
				    VirtualMatrix<T>& mat_interac,
				    const GlobalGenericMatrix<T>& nat_mat);
        
  public :
    void InitIndices(int n);
    void SetIndices(int i, const VectString& parameters);
    void SetPhysicalIndex(const string& name_media, int i, const VectString& parameters);
    Real_wp GetVelocityOfMedia(int ref) const;
    
    void ComputeRightHandSide(Vector<Complexe>& bsrc);
    
    void RunAll(const string& input_file, const string& name_elt);
    
    template<class T>
    void AddBoundaryTerms(VirtualMatrix<T>& mat_sp, GlobalGenericMatrix<T>& nat_mat);
    
    void ReadImpedance(int&, const VectString& param, int side);
    
  };
  
  
  //! 1-D gaussian source
  template<>
  class GaussianSource<Dimension1>
  {
  protected :
    Real_wp x0, radius;
    Real_wp alpha, beta;
    
  public :
    GaussianSource();

    Real_wp GetAmplitude(const Real_wp& x) const;
    
  };
  
  
  template<class T>
  class GaussianSourceField<T, Dimension1> : public GaussianSource<Dimension1>,
					     public VirtualSourceField<T, Dimension1>
  {
  public:
    void EvaluateFunction(const Real_wp& x, T& f) const;
    void EvaluateFunctionGradient(const Real_wp& x, T& f, TinyVector<T, 1>& grad_f) const;
    
  };
  

  //! class to solve Helmholtz equation in 1-D
  template<>
  class EllipticProblem<HelmholtzEquation1D> : public VarHelmholtz_1D<Complex_wp>
  {
  public :
    
  };


  //! class to solve Laplace equation in 1-D
  template<>
  class EllipticProblem<LaplaceEquation1D> : public VarHelmholtz_1D<Real_wp>
  {
  public :
    void ComputeElementaryMatrix(int i, IVect& num, VirtualMatrix<Real_wp>& A,
				 CondensationBlockSolver_Base<Real_wp>&,
				 const GlobalGenericMatrix<Real_wp>& nat_mat);
    
    void ComputeElementaryMatrix(int i, IVect& num, VirtualMatrix<Complex_wp>& A,
				 CondensationBlockSolver_Base<Complex_wp>&,
				 const GlobalGenericMatrix<Complex_wp>& nat_mat);
    
  };

} // namespace Montjoie

#define MONTJOIE_FILE_HELMHOLTZ_1D_HXX
#endif
