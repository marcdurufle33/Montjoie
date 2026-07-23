#ifndef MONTJOIE_FILE_TRANSMISSION_MODEL_ELASTIC_HXX

namespace Montjoie
{
  
  //! base class for transmission conditions and elastodynamic equation (see Burel's thesis)
  template<class Dimension>
  class VarTransmission_Elas_Base : public VarTransmission_Base<Dimension>
  {
  protected :
    bool curvilinear_case;

  public:
    Real_wp coef_ch, coef_sh;
    
  private:
    VarProblem<Dimension>& var_problem;
    VarBoundaryCondition_Base& var_boundary;

  public:
    template<class TypeEquation>
    VarTransmission_Elas_Base(EllipticProblem<TypeEquation>&);
    
    void InitTransmission();
    
    template<class Complexe>
    void AddTransmissionTermGen(const Complexe&, const Complexe&,
				const Complexe&, const Complexe&, const Complexe& alpha, 
				const GlobalGenericMatrix<Complexe>& nat_mat,
				VirtualMatrix<Complexe>& mat_sp, int, int);
    
  };


  template<class Complexe, class Dimension>
  class VarTransmission_Elas : public VarTransmission_Elas_Base<Dimension>
  {
  protected:
    //! physical properties of the thin layer
    Complexe rho_layer, lambda_layer, mu_layer, sigma_layer;
    
  private:
    VarProblem<Dimension>& var_problem;
    VarBoundaryCondition_Dim<Complexe, Dimension>& var_boundary;

  public:
    template<class TypeEquation>
    VarTransmission_Elas(EllipticProblem<TypeEquation>&);
    
    void SetInputData(const string&, const VectString& param);
    
    void AddTransmissionTerms(const Real_wp& alpha, 
			      const GlobalGenericMatrix<Real_wp>& nat_mat,
			      VirtualMatrix<Real_wp>& mat_sp, int, int);

    void AddTransmissionTerms(const Complex_wp& alpha, 
			      const GlobalGenericMatrix<Complex_wp>& nat_mat,
			      VirtualMatrix<Complex_wp>& mat_sp, int, int);

  };
  

  //! class for transmission conditions and elastodynamics (Burel's thesis)
  template<class Dimension>
  class VarTransmission<HarmonicElasticEquation<Dimension> > 
    : public VarTransmission_Elas<Complex_wp, Dimension>
  {
  public:
    template<class TypeEquation>
    VarTransmission(EllipticProblem<TypeEquation>&);

  };


  //! class for transmission conditions and elastodynamics (Burel's thesis)
  template<class Dimension>
  class VarTransmission<ElasticEquation<Dimension> > 
    : public VarTransmission_Elas<Real_wp, Dimension>
  {
  public:
    template<class TypeEquation>
    VarTransmission(EllipticProblem<TypeEquation>&);

  };
  
  
  //! base class for impedance operator in elastic transmission conditions
  /*!
    In 2-D
    \term \int coef_u u \cdot phi + coef_u_dphi u \cdot d_1 phi 
              + coef_phi_du phi \cdot d_1 u + coef_du_dphi d_1 u \cdot d_1 phi dx
    In 3-D
    term \int coef_u u \cdot phi + coef_u_d1_phi u \cdot d_1 phi + coef_u_d2_phi u \cdot d_2 phi dx
        + \int coef_phi_d1_u phi \cdot d1 u + coef_phi_d2_u phi \cdot d2 u dx
	+ \int coef_d1u_d1phi d1 u \cdot d1 phi + \int coef_d2u_d2phi d2 u \cdot d2 phi dx
	+ \int coef_d1u_d2phi d1 u \cdot d2 phi + \int coef_d1phi_d2u d1 phi \cdot d2 u dx
	
    u being the unknown (phi_j), and phi the test-function (phi_i)
   */
  template<class Dimension>
  class ImpedanceTransmission_Elas_Base
  {
  public :
    //! for the curvilinear case
    bool curvilinear, cond_on_gamma_plus;    
    
    //! normales and tangent vectors
    Vector<TinyVector<Real_wp, Dimension::dim_N> > VecNormale, VecTau, VecTauX, VecTauY;
    
  private:
    VarProblem<Dimension>& var_problem;
    VarTransmission_Elas_Base<Dimension>& var_transmission;
    
  public:
    ImpedanceTransmission_Elas_Base(VarProblem<Dimension>& var,
				    VarTransmission_Elas_Base<Dimension>& var_t);

    void SetFrenetBasis(const typename Dimension::R_N&, int);
    
    void ConstructFrenetBasis(int i, int num_elem, int num_edge, int num_loc, int k,
			      const SetPoints<Dimension>& Pts, const SetMatrices<Dimension>& Mat);
    
    template<class Complexe>
    void ApplyImpedancePhiGen(int n, int j, const TinyVector<Real_wp, 1>& phi_boundary,
			      const typename Dimension::R_N& dphi, Vector<Complexe>& feval_phi);

    template<class Complexe>
    void ApplyImpedanceGradGen(int n, int j, const TinyVector<Real_wp, 1>& phi_boundary,
			       const typename Dimension::R_N& dphi, Vector<Complexe>& feval_phi);
    
  };
  
  
  template<class Complexe, class Dimension>
  class ImpedanceTransmission_Elas : public ImpedanceTransmission_Elas_Base<Dimension>,
				     public ImpedanceFunction_Base<Complexe, Dimension>
  {
  public:
    //! coefficients mainly used in 2-D
    TinyMatrix<Complexe, General, Dimension::dim_N, Dimension::dim_N>
    coef_u, coef_du_dphi, coef_u_dphi, coef_phi_du;

    //! coefficients specific to the 3-D case
    TinyMatrix<Complexe, General, 3, 3> coef_u_d1_phi, coef_u_d2_phi, coef_phi_d1_u, coef_phi_d2_u,
      coef_d1u_d1phi, coef_d2u_d2phi, coef_d1u_d2phi, coef_d1phi_d2u;
    
  public :
    ImpedanceTransmission_Elas(VarProblem<Dimension>& var,
			       VarTransmission_Elas_Base<Dimension>& var_t);

    bool PresenceGradient() const;
    
    inline void EvaluateImpedancePhi(int i, int num_elem, int num_edge, int num_loc, int k,
				     const GlobalGenericMatrix<Complexe>& nat_mat, int ref,
				     const SetPoints<Dimension>& Pts,
				     const SetMatrices<Dimension>& Mat){}
    
    void EvaluateImpedanceGrad(int i, int num_elem, int num_edge, int num_loc, int k,
			       const GlobalGenericMatrix<Complexe>& nat_mat, int ref,
			       const SetPoints<Dimension>& Pts,
			       const SetMatrices<Dimension>& Mat);
    
    void ApplyImpedancePhi_H1(int n, int j, int, const TinyVector<Real_wp, 1>& phi_boundary,
                              const typename Dimension::R_N& dphi, Vector<Complexe>& feval_phi);

    void ApplyImpedanceGrad(int n, int j, int, const TinyVector<Real_wp, 1>& phi_boundary,
			    const typename Dimension::R_N& dphi, Vector<Complexe>& feval_phi);
 
  };

}

#define MONTJOIE_FILE_TRANSMISSION_MODEL_ELASTIC_HXX
#endif
