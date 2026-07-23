#ifndef MONTJOIE_FILE_IMPEDANCE_HELMHOLTZ_HXX

namespace Montjoie
{
  
  template<class Dimension>
  class VarGeneralizedImpedance_Helm;

  //! base class for first-order absorbing boundary condition for Helmholtz equation
  class ImpedanceABC_HelmholtzBase
  {
  public :
    //! type of absorbing boundary conditions
    enum{CLA_1=1, CLA_2=2, KM2_DELTA=3, KM2_OMEGA=4, KM2_OMEGA2=5, KM2_DELTA2=6,
         PARAM_KM1_DELTA=11, PARAM_KM1_OMEGA=12, PARAM_KM0_OMEGA=13,
         ATMO_1=21, ATMO_2=22, ATMO_3=23,
         ATMO_RBC_1=31, ATMO_SAI_0=32, ATMO_SAI_1=33, ATMO_HF_1=34, ATMO_ML_SPHERE=35, ATMO_ML_ANY=36,
         WHIT_1=37, WHIT_2=38, WHIT_3=39,
	 VISCOTHERMAL = 100, CREMER = 101};
    
    static void ComputeCoefficientsABC(Complex_wp myk, Real_wp courbure, int order, bool grazing_abc,
				       Complex_wp & alpha_cla, Complex_wp & beta_cla,
				       Complex_wp & delta_cla, Complex_wp& eta_cla,
				       Real_wp gamma_cla_coef, Real_wp theta_cla_coef, Real_wp zeta_cla_coef,
				       Real_wp h = 0.0, Real_wp g = 0.0);
    
    static int GetBoundaryConditionId(const IVect& ref, int pos, const VectString& parameters,
                                      bool& periodic, int& order, Real_wp& gamma, bool& take_curve);
    
  };
  
  
  //! class for first-order absorbing boundary condition for Helmholtz equation
  template<class Dimension>
  class ImpedanceABC_Helm : public ImpedanceFunction_Base<Complex_wp, Dimension>,
    public ImpedanceABC_HelmholtzBase
  {
  public :
    static bool gibc_condition;
    int num_row_phi, num_col_u;
    bool add_gibc_term;
    Complex_wp coef_abc, coef_abc_grad;
    Real_wp epsilon_pade, k_infty;
      
    const VarHelmholtzIndex_Base<Complex_wp, Dimension>& var_helm;  	
    const VarBoundaryCondition<Complex_wp, Dimension>& var_boundary;
    const VarProblem<Dimension>& var_problem;
    
  private :
    const VarGeneralizedImpedance_Helm<Dimension>& var_gibc;
    
    void InitDefaultValues();
  public :
      
    template<class TypeEquation>
    ImpedanceABC_Helm(const EllipticProblem<TypeEquation>& var);

    bool PresenceGradient() const;
    
    Complex_wp GetCoefficient(int i, int iquad, int num_loc, int k, int ref_d, int ref,
			      const SetPoints<Dimension>&, const SetMatrices<Dimension>&);
    
    // computation of impedance
    void EvaluateImpedancePhi(int i, int iquad, int num_edge, int num_loc, int k,
			      const GlobalGenericMatrix<Complex_wp>& nat_mat, int ref_domain,
			      const SetPoints<Dimension>&, const SetMatrices<Dimension>&);
    
    void EvaluateImpedanceGrad(int i, int iquad, int num_edge, int num_loc, int k,
			       const GlobalGenericMatrix<Complex_wp>& nat_mat, int ref_domain,
			       const SetPoints<Dimension>&, const SetMatrices<Dimension>&);
        
  };


  //! class for first-order absorbing boundary condition for Helmholtz equation
  template<class Dimension>
  class ImpedanceABC<Complex_wp, HelmholtzEquation<Dimension> > 
    : public ImpedanceABC_Helm<Dimension>
  {    
  public :
    ImpedanceABC(const EllipticProblem<HelmholtzEquation<Dimension> >& var);
    
  };


  //! class for first-order absorbing boundary condition for Helmholtz equation
  template<class Dimension>
  class ImpedanceABC<Complex_wp, HelmholtzEquationDG<Dimension> > 
    : public ImpedanceABC_Helm<Dimension>
  {    
  public :
    ImpedanceABC(const EllipticProblem<HelmholtzEquationDG<Dimension> >& var);
    
  };


  template<class Dimension>
  bool ImpedanceABC_Helm<Dimension>:: gibc_condition(false);
  
#ifdef MONTJOIE_WITH_HIGH_CONDUCTIVITY_MODEL
  //! variable impedance for high-conductivity condition
  template<class Dimension>
  class ImpedanceHighConductivity<Complex_wp, HelmholtzEquation<Dimension> >
    : public ImpedanceFunction_Base<Complex_wp, Dimension>
  {
  private :
    const VarHelmholtzIndex_Base<Complex_wp, Dimension>& var_helm;
    const VarBoundaryCondition_Dim<Complex_wp, Dimension>& var_boundary;
    const VarProblem<Dimension>& var_problem;
    
  public :
    static Complex_wp alpha; //!< sqrt(i)

    ImpedanceHighConductivity(const EllipticProblem<HelmholtzEquation<Dimension> >& var);

    bool PresenceGradient() const;
    
    void EvaluateImpedancePhi(int i, int num_elem, int num_edge, int num_loc, int k,
			      const GlobalGenericMatrix<Complex_wp>& nat_mat, int ref_d,
			      const SetPoints<Dimension>& Pts,
			      const SetMatrices<Dimension>& MatricesElem);
    
    void EvaluateImpedanceGrad(int i, int num_elem, int num_edge, int num_loc, int k,
			       const GlobalGenericMatrix<Complex_wp>& nat_mat, int ref_d,
			       const SetPoints<Dimension>& Pts, const SetMatrices<Dimension>& Mat);
    
  };

  
#endif  


  //! base class for GIBC and Helmholtz equation
  template<class Dimension>
  class VarGeneralizedImpedance_Helm : public VarGeneralizedImpedance_Base
  {
  private :
    typedef Vector<typename Dimension::R_N> VectR_N;
    VarProblem<Dimension>& var_problem;
    VarHelmholtzIndex_Base<Complex_wp, Dimension>& var_helm;
    VarBoundaryCondition<Complex_wp, Dimension>& var_boundary;
    
  public:
    template<class TypeEquation>
    VarGeneralizedImpedance_Helm(EllipticProblem<TypeEquation>& var) ;
    
    void TreatGibc(const IVect& Epart);
    
    void InitGIBC();

    void AddGibcTerms(const Complex_wp& alpha, const GlobalGenericMatrix<Complex_wp>& nat_mat,
                      VirtualMatrix<Complex_wp>& mat_sp, int offset_row = 0, int offset_col = 0);

    void AddGibcTerms(const Real_wp& alpha, const GlobalGenericMatrix<Real_wp>& nat_mat,
                      VirtualMatrix<Real_wp>& mat_sp, int offset_row = 0, int offset_col = 0);
    
  };


  //! class for GIBC and Helmholtz equation
  template<class Dimension>
  class VarGeneralizedImpedance<HelmholtzEquation<Dimension> > 
    : public VarGeneralizedImpedance_Helm<Dimension>
  {    
  public :
    VarGeneralizedImpedance(EllipticProblem<HelmholtzEquation<Dimension> >& var);
    
  };


  template<class Dimension>
  class VarGeneralizedImpedance<HelmholtzEquationDG<Dimension> > 
    : public VarGeneralizedImpedance_Helm<Dimension>
  {    
  public :
    VarGeneralizedImpedance(EllipticProblem<HelmholtzEquationDG<Dimension> >& var);
    
  };

}

#define MONTJOIE_FILE_IMPEDANCE_HELMHOLTZ_HXX
#endif

