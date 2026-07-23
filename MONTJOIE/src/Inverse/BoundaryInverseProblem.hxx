#ifndef MONTJOIE_FILE_BOUNDARY_INVERSE_PROBLEM_HXX

namespace Montjoie
{

  template<class Dimension>
  class BoundaryInverseProblemMultiFreqParam;
  
  //! Base class for boundary inverse problems
  template<class Dimension>
  class BoundaryInverseProblem_Base : public itreg::ForwardOperator<Real_wp>,
                                      public MultiFrequencyProblem
  {
    template<class Dim>
    friend class BoundaryInverseProblemMultiFreqParam;

  public:
    typedef typename Dimension::R_N R_N;
    typedef Real_wp (*type_ptr_fctC)(const R_N&);

  protected:
    typedef typename Dimension::VectR_N VectR_N;
    

    //! functions defining the exact parameters c
    Vector<type_ptr_fctC> c_exact;
  
    //! object for solving the time-harmonic equation
    VarHarmonic_Base<Complex_wp, Dimension>& var_helm;
    
    //! object solving the linear systems
    Vector<All_LinearSolver*> linear_solver;
    
    //! stored solutions (for all explored frequencies)
    Vector<Matrix<Complex_wp, General, ColMajor> > sol_u;
    Vector<Matrix<Complex_wp, General, ColMajor> > u_observable;
    
    //! stored sources
    Matrix<Complex_wp, General, ColMajor> source_rhs;
    
    //! exact values of all parameters C stored
    VectReal_wp C_store_exact, C_init_guess;
  
    //! object discretizing parameters c
    Vector<VariableParameterND_Base<Real_wp, Dimension>* > c_discrete;
    Vector<string> param_c_name; Vector<int> param_c_type;
    
    // parameters for the source
    int nb_comp_source_invert;
    int nb_angles_measure, nb_angles_to_invert;
    VectReal_wp coef_combination_measure;
    
    // quadrature points where parameters c needs to be evaluated
    Vector<Vector<R_N> > points_quad;
    Vector<VectReal_wp> eval_c_quad;
    
    //! matrix interpolating c on quadrature points
    Matrix<Real_wp, General, ArrayRowSparse> CollMatrix;
  
    //! parameters associated with the regularization method
    Vector<string> regularization_parameter;
    Vector<string> approx_c_parameter;
    
    // Dirichlet condition
    Vector<Vector<string> > dirichlet_c;
    Vector<IVect> DofDirichletC;
    Vector<VectReal_wp> ValueDirichletC;
    Vector<bool> elt_to_keep;
    
    VectReal_wp coef_input_frequency;
    
    //! object defining sources and receptors
    VarMigration_Base<Dimension>& var_migration;
    
    //! reference of the elements where c is discretized
    IVect ref_volume;
    
    //! weights for measurements
    VectReal_wp weights_measure;
    int measured_part, nb_observables;
    enum {ALL_PART, REAL_PART, IMAG_PART, MODULUS_PART};
    
    //! weights for basis functions
    Vector<VectReal_wp> weights_basis;
    VectReal_wp param_scaling;
    
    // finite element used and equation
    string name_element, name_equation;

    // experimental data
    bool write_experimental_data;
    string file_name_experimental_data;    
    int NdisplayC;

    int number_facto_lu, number_solve;
    
  public:

    //! array of static functions
    static type_ptr_fctC exact_fct_to_find[42];
    
    BoundaryInverseProblem_Base(VarHarmonic_Base<Complex_wp, Dimension>& var,
				VarMigration_Base<Dimension>& var_mig);

    ~BoundaryInverseProblem_Base();
    
    inline VarMigration_Base<Dimension>& GetMigrationProblem() { return var_migration; }
    
    void ResizeNbParameters(int n);
    virtual void SetExactParameter(int k, const string& name, type_ptr_fctC fct);
    
    const Vector<string>& GetRegularizationParameter() const;
    void SetInputData(const string& keyword, const Vector<string>& param);

    bool ExactSolutionKnown() const;
    Real_wp DistanceToExactSolution(const VectReal_wp& x);
    
    void ExtractParametersC(const VectReal_wp& cn0, Vector<VectReal_wp>& cn);
    void CopyParametersC(const VectReal_wp& cn0, Vector<VectReal_wp>& cn);

    Real_wp GetInterpolationError();    
    void WriteDatas(const VectReal_wp& cn0);
    
    void Construct(const Vector<string>& lines_montjoie, const string& name_elt,
		   const string name_eq);

    void ProvideExactValues(const VectReal_wp& exact_cte);
    
    void ConstructDirichlet();
    void ImposeDirichletCondition(Vector<VectReal_wp>& cn);

    void FindInitGuess(VectReal_wp& cn_init);
    void SetVariableIndex(const Vector<VectR_N>&,
			  const Vector<VectReal_wp>& eval_c,
			  VarHarmonic_Base<Complex_wp, Dimension>& var);

    void FillVariableC(const VectReal_wp& cn0);

    void ModifySource(const VarHarmonic_Base<Complex_wp, Dimension>& var,
                      const VectReal_wp& cn0,
                      Matrix<Complex_wp, General, ColMajor>& source);
    
    void InitCombinationMeasure(const VectReal_wp& cn0);

    void Init(const VectReal_wp& cn0);
    void Evaluate(VectReal_wp& y);
    
    
    void MltParamDeriv(const Real_wp& alpha, const Vector<VectReal_wp>& c_quad,
		       Vector<Matrix<Complex_wp, General, ColMajor> >& y);

    void MltParamDerivOmega(const Real_wp& alpha, int k, const Vector<VectReal_wp>& c_quad,
			    const Matrix<Complex_wp, General, ColMajor>& x,
			    Matrix<Complex_wp, General, ColMajor>& y);
    
    void MltParamDerivTrans(const Real_wp& alpha, const Vector<Matrix<Complex_wp, General, ColMajor> >& y,
			    Vector<VectReal_wp>& ProdC);
  
    void Derivative(const VectReal_wp& h0, VectReal_wp& res);
    void Adjoint(const VectReal_wp& G, VectReal_wp& res);

    void EvaluateJacobian(const VectReal_wp& x, VectReal_wp& f, Matrix<Real_wp>& df);
    
    void EvaluateDifferently(VectReal_wp& y);
    
    void AppendData(const Matrix<Complex_wp, General, ColMajor>& data,
		    int& num, VectReal_wp& y);

    virtual void CombineData(const Matrix<Complex_wp, General, ColMajor>& data0,
                             Matrix<Complex_wp, General, ColMajor>& data);
    
    virtual void CombineDeriveData(const VectReal_wp& h0, int k,
                                   const Matrix<Complex_wp, General, ColMajor>& data_simu,
                                   Matrix<Complex_wp, General, ColMajor>& data_comb);

    virtual void CombineAdjointData(int k, const Matrix<Complex_wp, General, ColMajor>& data_simu,
                                    Matrix<Complex_wp, General, ColMajor>& data_comb, VectReal_wp& h0);
    
    void ScatterData(int k, const Matrix<Complex_wp, General, ColMajor>& data_simu,
                     int& num, VectReal_wp& res);
    
    void ExtractData(int k, const VectReal_wp& y,
		     int& num, Matrix<Complex_wp, General, ColMajor>& data);

    void WriteMeasuredY(const VectReal_wp& y, const string& output_name);
    
    void GiveIterate(int n, const VectReal_wp& x);

    virtual void InitIndices() = 0;
    
    virtual void SetInversionParameters(const VectString& param) = 0;
    
    virtual void GetVariableIndex(int ref, VarHarmonic_Base<Complex_wp, Dimension>& var,
				  Vector<PhysicalVaryingMedia<Dimension, Real_wp>* >&,
				  Vector<PhysicalVaryingMedia<Dimension, Complex_wp>* >&) = 0;

    virtual void SetVariableIndexPoint(int ref, int i, int i1, int j, int offset, const Vector<VectReal_wp>& eval_c,
				       Vector<PhysicalVaryingMedia<Dimension, Real_wp>* >& rho_real,
				       Vector<PhysicalVaryingMedia<Dimension, Complex_wp>* >& rho_cplx) = 0;
    
    virtual void MltParamDerivLocal(int k, const Real_wp& omega2, VectComplex_wp& u_glob,
				    int i, const Real_wp& alpha, const Vector<VectReal_wp>& c_quad,
				    int offset, VectComplex_wp& y_col) = 0;

    virtual void MltParamDerivTransLoc(int k, const Real_wp& omega2, VectComplex_wp& u_glob,
				       int i, const Real_wp& alpha, int offset, 
				       const VectComplex_wp& y_col, Vector<VectReal_wp>& ProdC) = 0;
    
  };


  //! class for solving multi-frequency for EvaluateDifferently
  template<class Dimension>
  class BoundaryInverseProblemMultiFreqParam : public MultiFrequencyParameter
  {
    template<class Dim> friend class BoundaryInverseProblem_Base;
    
  protected:
    //! object defining sources and receptors
    VarMigration_Base<Dimension>& var_migration;
    // boundary inverse problem
    BoundaryInverseProblem_Base<Dimension>& var_inverse;
    // synthetic data
    VectReal_wp& y;

    int num;
    Matrix<Complex_wp, General, ColMajor> source_rhs_ref;
    Matrix<Complex_wp, General, ColMajor> data_exp;
    
  public:
    BoundaryInverseProblemMultiFreqParam(BoundaryInverseProblem_Base<Dimension>& var, VectReal_wp&);
    
    void SetSourceFrequency(int k, int kg, int m, int mg,
                            const VarComputationProblem&, Matrix<Complex_wp, General, ColMajor>& sol);
    void GiveSolutionFrequency(int k, int kg, int m, int mg,
                               const VarComputationProblem&, Matrix<Complex_wp, General, ColMajor>& sol_exp);

  };
  
  
  template<class TypeEquation>
  class BoundaryInverseProblem : public BoundaryInverseProblem_Base<typename TypeEquation::Dimension>
  {
  };
  
  
  //! Inverse problems for 2-D or 3-D Helmholtz equation
  template<class Dimension>
  class BoundaryInverseProblem<HelmholtzEquation<Dimension> > : public BoundaryInverseProblem_Base<Dimension>
  {
  protected:
    EllipticProblem<HelmholtzEquation<Dimension> > var_helm;
    VarMigration<HelmholtzEquation<Dimension> > var_migration;
    
    typedef typename Dimension::R_N R_N;
    typedef Real_wp (*type_ptr_fctC)(const R_N&);
    
    // available parameters (first parameters are NATURAL, last ones are HELIO)
    enum {PARAM_RHO, PARAM_MU, PARAM_MX, PARAM_MY, PARAM_MZ, PARAM_SIGMA,
          PARAM_DENSITY, PARAM_C, PARAM_UX, PARAM_UY, PARAM_UZ, PARAM_GAMMA};
    
    int type_convention;
    enum {NONE, NATURAL, HELIO};
    
    bool presence_mass, presence_stiff, presence_flow;
    
    VectReal_wp ref_density, ref_c, ref_gamma;
    
    Vector<VectReal_wp> variable_density, variable_c, variable_gamma;
    Vector<typename Dimension::VectR_N> variable_M;

    VectReal_wp ref_rho, ref_mu, ref_sigma, ref_beta;
    typename Dimension::VectR_N ref_M;

    Vector<VectReal_wp> variable_rho, variable_mu, variable_sigma, variable_beta;
    
  public:
    BoundaryInverseProblem();

    void InitIndices();

    void SetInversionParameters(const VectString& param);    
    void SetExactParameter(int k, const string& name, type_ptr_fctC fct);
    int GetIdNumberParameter(const string& name) const;
    
    void GetVariableIndex(int ref, VarHarmonic_Base<Complex_wp, Dimension>& var,
			  Vector<PhysicalVaryingMedia<Dimension, Real_wp>* >&,
			  Vector<PhysicalVaryingMedia<Dimension, Complex_wp>* >&);

    void SetVariableIndexPoint(int ref, int i, int i1, int j, int offset, const Vector<VectReal_wp>& eval_c,
			       Vector<PhysicalVaryingMedia<Dimension, Real_wp>* >& rho_real,
			       Vector<PhysicalVaryingMedia<Dimension, Complex_wp>* >& rho_cplx);

    void GetPhysicalIndexesElem(int i, int offset, int nb_pts, TinyVector<VectReal_wp, Dimension::dim_N>& flow,
                                VectReal_wp& beta, VectReal_wp& rho, VectReal_wp& c, VectReal_wp& gamma);
    
    void MltParamDerivLocal(int k, const Real_wp& omega2, VectComplex_wp& u_glob,
			    int i, const Real_wp& alpha, const Vector<VectReal_wp>& c_quad,
			    int offset, VectComplex_wp& y_col);
    
    void MltParamDerivTransLoc(int k, const Real_wp& omega2, VectComplex_wp& u_glob,
			       int i, const Real_wp& alpha, int offset, 
			       const VectComplex_wp& y_col, Vector<VectReal_wp>& ProdC);
    
  };


  //! Inverse problems for axisymmetric Helmholtz equation
  template<>
  class BoundaryInverseProblem<HelmholtzEquationAxi> : public BoundaryInverseProblem_Base<Dimension2>
  {
  protected:
    EllipticProblem<HelmholtzEquationAxi> var_helm;
    VarMigration<HelmholtzEquationAxi> var_migration;
    
    typedef Real_wp (*type_ptr_fctC)(const R2&);
    
    // available parameters
    enum {PARAM_RHO, PARAM_MU, PARAM_DENSITY, PARAM_C};
    
    VectComplex_wp ref_density, ref_c, ref_rho, ref_mu, ref_sigma;
    
  public:
    BoundaryInverseProblem();

    void InitIndices();

    void SetInversionParameters(const VectString& param);
    void SetExactParameter(int k, const string& name, type_ptr_fctC fct);
    
    void GetVariableIndex(int ref, VarHarmonic_Base<Complex_wp, Dimension2>& var,
			  Vector<PhysicalVaryingMedia<Dimension2, Real_wp>* >&,
			  Vector<PhysicalVaryingMedia<Dimension2, Complex_wp>* >&);

    void SetVariableIndexPoint(int ref, int i, int i1, int j, int offset, const Vector<VectReal_wp>& eval_c,
			       Vector<PhysicalVaryingMedia<Dimension2, Real_wp>* >& rho_real,
			       Vector<PhysicalVaryingMedia<Dimension2, Complex_wp>* >& rho_cplx);

    void MltParamDerivLocal(int k, const Real_wp& omega2, VectComplex_wp& u_glob,
			    int i, const Real_wp& alpha, const Vector<VectReal_wp>& c_quad,
			    int offset, VectComplex_wp& y_col);
    
    void MltParamDerivTransLoc(int k, const Real_wp& omega2, VectComplex_wp& u_glob,
			       int i, const Real_wp& alpha, int offset, 
			       const VectComplex_wp& y_col, Vector<VectReal_wp>& ProdC);
    
  };
  
  
  //! Inverse problems for 2-D or 3-D time-harmonic elastic equation
  template<class Dimension>
  class BoundaryInverseProblem<HarmonicElasticEquation<Dimension> > : public BoundaryInverseProblem_Base<Dimension>
  {
  protected:
    EllipticProblem<HarmonicElasticEquation<Dimension> > var_helm;
    VarMigration<HarmonicElasticEquation<Dimension> > var_migration;

    typedef typename Dimension::R_N R_N;
    typedef Real_wp (*type_ptr_fctC)(const R_N&);

    // available parameters
    enum {PARAM_RHO, PARAM_LAMBDA, PARAM_MU, PARAM_LAMBDA_I, PARAM_MU_I,
	  PARAM_C11, PARAM_C22, PARAM_C33, PARAM_C12, PARAM_C13,
	  PARAM_C23, PARAM_C44, PARAM_C55, PARAM_C66,
	  PARAM_C14, PARAM_C15, PARAM_C16, PARAM_C24, PARAM_C25, PARAM_C26,
	  PARAM_C34, PARAM_C35, PARAM_C36, PARAM_C45, PARAM_C46, PARAM_C56,
	  PARAM_C11_I, PARAM_C22_I, PARAM_C33_I, PARAM_C12_I, PARAM_C13_I, 
	  PARAM_C23_I, PARAM_C44_I, PARAM_C55_I, PARAM_C66_I,
	  PARAM_C14_I, PARAM_C15_I, PARAM_C16_I, PARAM_C24_I, PARAM_C25_I, PARAM_C26_I,
	  PARAM_C34_I, PARAM_C35_I, PARAM_C36_I, PARAM_C45_I, PARAM_C46_I, PARAM_C56_I};
    
    Vector<Complex_wp> ref_density;
    Vector<TinySymmetricTensor<Complex_wp, Dimension::dim_N> > ref_C;
    
  public:
    BoundaryInverseProblem();

    void InitIndices();

    void SetInversionParameters(const VectString& param);
    void SetExactParameter(int k, const string& name, type_ptr_fctC fct);

    void CombineData(const Matrix<Complex_wp, General, ColMajor>& data0,
                             Matrix<Complex_wp, General, ColMajor>& data);
    
    void CombineDeriveData(const VectReal_wp& h0, int k,
                           const Matrix<Complex_wp, General, ColMajor>& data_simu,
                           Matrix<Complex_wp, General, ColMajor>& data_comb);

    void CombineAdjointData(int k, const Matrix<Complex_wp, General, ColMajor>& data_simu,
                            Matrix<Complex_wp, General, ColMajor>& data_comb, VectReal_wp& h0);
    
    static void GetIndexParam(int num, int& i, int& j, int& k, int& l, bool& is_imag);
    
    void GetVariableIndex(int ref, VarHarmonic_Base<Complex_wp, Dimension>& var,
			  Vector<PhysicalVaryingMedia<Dimension, Real_wp>* >&,
			  Vector<PhysicalVaryingMedia<Dimension, Complex_wp>* >&);

    void SetVariableIndexPoint(int ref, int i, int i1, int j, int offset, const Vector<VectReal_wp>& eval_c,
			       Vector<PhysicalVaryingMedia<Dimension, Real_wp>* >& rho_real,
			       Vector<PhysicalVaryingMedia<Dimension, Complex_wp>* >& rho_cplx);
    
    void MltParamDerivLocal(int k, const Real_wp& omega2, VectComplex_wp& u_glob,
			    int i, const Real_wp& alpha, const Vector<VectReal_wp>& c_quad,
			    int offset, VectComplex_wp& y_col);
    
    void MltParamDerivTransLoc(int k, const Real_wp& omega2, VectComplex_wp& u_glob,
			       int i, const Real_wp& alpha, int offset, 
			       const VectComplex_wp& y_col, Vector<VectReal_wp>& ProdC);
    
  };
  
}

#define MONTJOIE_FILE_BOUNDARY_INVERSE_PROBLEM_HXX
#endif
