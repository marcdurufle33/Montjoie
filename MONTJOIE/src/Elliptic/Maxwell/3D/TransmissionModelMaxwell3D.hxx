#ifndef MONTJOIE_FILE_TRANSMISSION_MODEL_MAXWELL3D_HXX

namespace Montjoie
{

  class VarTransmission_Maxwell3D;
  
  //! base class for impedance boundary condition
  /*!
    if var_transmission.homogeneous_thin_layer is false
    we add the terms
    \int_\Gamma C1 n \times \varphi_i \cdot \varphi_j  dx
    - \int_\Gamma n \times (C2 n \times \varphi_i) \cdot \varphi_j  dx
    + \int_\Gamma C3 n \cdot curl \varphi_i   n \cdot curl \varphi_j  dx
   */
  class ImpedanceTransmissionMaxwell3D
    : public ImpedanceFunction_Base<Complex_wp, Dimension3>
  {
  public :
    Complex_wp C1, C3 ; 
    TinyMatrix<Complex_wp, Symmetric, 3, 3> C2 ; 
    Vector<TinyVector<Complex_wp, 3> > invMu_normale;
    //! coefficients used if var_transmission.homogeneous_thin_layer is true
    /*!
      we add the term
      coef_val \int curl(\varphi_i) \ times n/mu \cdot \varphi_j
      + coef_grad \int n/mu \times \varphi_i \cdot curl(\varphi_j)
      i.e
      - coef_val \int n/mu \times curl(U) \cdot E
      - coef_grad \int n/mu \times curl(E) \cdot U
     */
    Complex_wp coef_val, coef_grad, coef_u;

  private:
    Vector<TinyMatrix<Complex_wp, General, 3, 3> > feval, feval_grad;
    const HarmonicMaxwell_3D<Complex_wp>& var_maxwell;
    const VarTransmission_Maxwell3D& var_transmission;
    
  public:
    template<class TypeEquation>
    ImpedanceTransmissionMaxwell3D(const EllipticProblem<TypeEquation>& vars);
    
    void EvaluateImpedancePhi(int i, int iquad, int num_edge, int num_loc, int k,
			      const GlobalGenericMatrix<Complex_wp>& nat_mat, int ref_domain,
			      const SetPoints<Dimension3>&, const SetMatrices<Dimension3>&);
    
    void EvaluateImpedanceGrad(int i, int iquad, int num_edge, int num_loc, int k,
			       const GlobalGenericMatrix<Complex_wp>& nat_mat, int ref_domain,
			       const SetPoints<Dimension3>&, const SetMatrices<Dimension3>&);

    bool PresenceGradient() const;
    
    // computes f(phi_i, dphi_i) (m = unknown number associated with row i)
    void ApplyImpedancePhi_Hcurl(int m, int j, int, const TinyVector<Real_wp, 3>& phi,
                                 const TinyVector<Real_wp, 3>& curl_phi, VectComplex_wp& f_phi);
    
    // computes g(phi_i, dphi_i) (m = unknown number associated with row i)
    void ApplyImpedanceCurl(int m, int j, int, const TinyVector<Real_wp, 3>& phi,
			    const TinyVector<Real_wp, 3>& curl_phi, VectComplex_wp& g_phi);
    
  };  

  
  //! base class for transmission conditions for 3-D Maxwell's equations
  /*!
    The transmission conditions are written by introducing an additional unknown lambda :
    \lambda = < 1/mu curl E >_alpha
    where <u>_alpha = alpha u^+ + (1-alpha) u^-
    [u] = u^+ - u^-
    then the transmission conditions are given as
    [n \times E] = delta (-A_\beta curl curl \lambda + B_\beta \lambda)
    [n/mu \times curl E] = delta (-C_\beta curl curl <E_T>_{1-alpha} + D_\beta <E_T>_{1-alpha})
    where delta is the thickness of the layer
    A_\beta is stored in coefficient A, B_beta in coefficient D1 (which can be a matrix),
    C_\beta in coefficient B and D_\beta in coefficient D2 (which can be a matrix)
    
    The variational formulation of these models is equal to
    -omega^2 \int_\Omega \epsilon E \cdot U dx + \int_\Omega 1/mu curl E \cdot curl U dx
    - \int_\Gamma n \times \lambda \cdot [U] dx
    + \delta \int_\Gamma C_\beta curl <E_T>_{1-\alpha} curl <U_T>_{1-\alpha} dx
    - \delta \int_\Gamma D_\beta <E_T>_{1-\alpha} \cdot <U_T>_{1-\alpha} dx = source
    
    \int_\Gamma [n \times E] \cdot \xi dx 
    + \delta \int_\Gamma A_\beta curl \lambda curl \xi dx
    - delta \int_Gamma B_beta \lambda \cdot \xi = 0
    
    Such variational formulation is used for the following models :
    MODEL_DELOURME : uniform thin layer, Gamma is the median surface (model denoted UNIFORM_SYM in the data file)
    MODEL_PERON : uniform thin layer, Gamma is the surface \Gamma- (model denoted UNIFORM_UNSYM in the data file)
    MODEL_PERON_RESISTIVE : uniform thin layer resistive (high constrast between the conductivity
    outside and inside the layer, in which the conductivity drops, model denoted RESISTIVE)
    MODEL_PERON_RESISTIVE_SIGMA : uniform thin layer resistive (high constrast between the conductivity
    outside and inside the layer, in which the conductivity drops, model denoted RESISTIVE_SIGMA)
    MODEL_CHUN : uniform thin layer, Gamma is thick, Decalage should be set to 0.5
    MODEL_PARAMETRIC : periodic thin layer (see Delourme's thesis), the coefficients A, B, D1, D2
    are directly given in the data file. They can be pre-computed by calling CellulairePeriodique3D
    
    For the conductive model (MODEL_CONDUCTIVE), two additional unknowns are considered
    \lambda_1 = < 1/mu curl E >   (\alpha = 1/2)
    \lambda_2 = [ 1/mu curl E ]
    Transmission conditions are given as
    [E \times n] = \delta (L_1 \lambda_1 + L_3 \lambda_2)
    <E \times n> = \delta (L_3 \lambda_1 + L_2 \lambda_2)
    with L_i = A_i curl curl - B_i
    coefficients A1, A2, A3, B1, B2, B3 are equal to :
    A1 = -1/omega^2 < (epsilon + i sigma/omega)^{-1} >
    B1 = 2 mu_0 / gamma tanh(gamma/2) - <mu>
    A2 = A1 / 4
    B2 = mu_0/(2 gamma) cotanh(gamma/2) - <mu> / 4
    A3 = - 1/(4 omega^2) [ (epsilon + i sigma/omega)^{-1} ]
    B3 = - 1/4 [mu]
    with gamma = exp(3i pi / 4)
   */
  class VarTransmission_Maxwell3D : public VarTransmission_Base<Dimension3>
  {
    friend class ImpedanceTransmissionMaxwell3D;

  protected :
    // parameters of the thin layer
    // penalization is present to enforce unicity of the linear system
    Complex_wp epsilon_layer, mu_layer, sigma_layer, penalization_layer;
    bool homogeneous_thin_layer; // if true, a test case is considered
    bool resistive_model; // no longer used
    Real_wp alpha; // parameter in the variational formulation (alpha=0.5 -> centered)
    Complex_wp A, B ; // parameters of the transmission conditions
    int type_model_thin_layer; // model considered
    // available models :
    enum {MODEL_DELOURME, MODEL_PERON, MODEL_PERON_RESISTIVE,
          MODEL_PERON_RESISTIVE_SIGMA, MODEL_PARAMETRIC, MODEL_CHUN, MODEL_CONDUCTIVE};
    
    // parameters of the transmission conditions
    TinyMatrix<Complex_wp, Symmetric, 3, 3> D1, D2;
    
    // parameters for a conductive material
    Complex_wp A1, B1, A2, B2, A3, B3;
    
    ImpedanceTransmissionMaxwell3D fimpedance;

  private:
    VarProblem<Dimension3>& var_problem;
    HarmonicMaxwell_3D<Complex_wp>& var_maxwell;
    VarBoundaryCondition_Dim<Complex_wp, Dimension3>& var_boundary;
    
    void InitDefaultValues();

  public:
    template<class TypeEquation>
    VarTransmission_Maxwell3D(EllipticProblem<TypeEquation>&) ;

    bool IsSymmetricCondition() const;
    
    void SetInputData(const string&, const VectString& param);
    
    void UpdateNumberOfDofs(int& nodl, int& nodl_scalar);
    void UpdateGlobalDofs();
    
    void InitTransmission();

    void AddTransmissionTerms(const Real_wp& alpha,
			      const GlobalGenericMatrix<Real_wp>& nat_mat,
			      VirtualMatrix<Real_wp>& mat_sp, int, int);
    
    void AddTransmissionTerms(const Complex_wp& alpha,
			      const GlobalGenericMatrix<Complex_wp>& nat_mat,
			      VirtualMatrix<Complex_wp>& mat_sp, int, int);
    
  };
  
  
  //! class for transmission conditions for 3-D Maxwell's equations
  template<>
  class VarTransmission<HarmonicMaxwellEquation_3D> 
    : public VarTransmission_Maxwell3D
  {
  public:
    template<class TypeEquation>
    VarTransmission(EllipticProblem<TypeEquation>&);
    
  };
  
} // end namespace

#define MONTJOIE_FILE_TRANSMISSION_MODEL_MAXWELL3D_HXX
#endif

  
