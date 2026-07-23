#ifndef MONTJOIE_FILE_BOUNDARY_CONDITION_HARMONIC_HXX

namespace Montjoie
{
  template<class T, class Dimension>
  class ImpedanceFunction_Base
  {
  };

  
  //! base class for treatment of boundary conditions
  class VarBoundaryCondition_Base
  {
  public :           
    //! column numbers to use for impedance matrix
    IVect NewColumnNumbers_Impedance, ProcColumnNumbers_Impedance;
    //! row numbers to use for impedance matrix
    IVect NewRowNumbers_Impedance, ProcRowNumbers_Impedance;
    
  protected :
    
    /***********************
     * Dirichlet condition *
     ***********************/
    
    //! array of degrees of freedom numbers with a dirichlet condition 
    IVect Dirichlet_dof;
    Vector<bool> is_dof_dirichlet;
    //! number of dofs with Dirichlet condition
    int nb_dof_dirichlet, nb_dof_dirichlet_all;  bool matrix_symmetric_dir;
    Real_wp coef_dirichlet_matrix;
    
    //! which components are set to 0 for supported condition
    Vector<IVect> supported_components_BC;
    
    // order_high_conductivity is used for high conductivity objects
    // order 1 is the case of the perfect conductor
    // order 2, 3 and 4 are implemented
    int order_high_conductivity; //!< order of high conductivity boundary condition
    
    
    /**********************
     * Periodic condition *
     **********************/

    //! mode numbers to compute
    IVect list_number_mode;
    int number_mode; //!< current mode number
    TinyVector<int, 3> nb_modes_periodic_xyz;
    bool plane_wave_quasi_periodic;
    //! type of symmetry (in x, y, z, theta, etc)
    int type_symmetry;
    //! if true, the modes are not stored
    bool do_not_store_modes;
    
  public :
    enum {NO_SYMMETRY, PERIODIC_THETA, PERIODIC_X, PERIODIC_Y, PERIODIC_Z,
          PERIODIC_XY, PERIODIC_XZ, PERIODIC_YZ, PERIODIC_XYZ, PERIODIC_ZTHETA};

    
    /***************
     * ABC and PML *
     ***************/
    

    //! take into account curvature for absorbing condition ?
    bool take_into_account_curvature_for_abc, grazing_abc;
    // ! coef of the parametrized ABC - to be multiplied by the curvature
    Real_wp  gamma_cla_coef, theta_cla_coef, zeta_cla_coef;

  protected :    
    //! order of the absorbing boundary condition
    int order_ABC;
    
    int nb_eltPML, nb_eltPML_all; //!< number of elements inside PML
    // parameters of the PML (Perfectly Matched Layers)
    Real_wp vsigma; //!< damping factor for PML
    Real_wp max_velocity_pml; //!< maximum speed in PML
    
    // damping function (parabole, constant, etc)
    int function_damping_pml;
    enum{ PML_PARABOLE, PML_CONSTANT, PML_SHIFTED_PARABOLE, PML_LINEAR};
    Real_wp offset_damping_pml, coef_mixed_damping_pml;

    //Real_wp omega_pml;
    //bool dispersive_pml;
    
  private :
    VarProblem_Base& var_problem;
    VarComputationProblem& var_computation;
    DistributedProblem_Base& var_comm;
    VarSourceProblem_Base& var_source;
    
    void InitDefaultValues();
    
  public :
    
    template<class TypeEquation>
    VarBoundaryCondition_Base(EllipticProblem<TypeEquation>&);
    
    void SetInputData(const string& description_field, const VectString& parameters);
    virtual bool GetInitialSymmetrization() const;
    
    virtual int GetBoundaryConditionId(const IVect&, int pos, const VectString&, bool& per);
    
  protected :
    virtual void CheckAndReallocateReferences(int ref) = 0;
    virtual void FillParameterCondition(const Vector<string>& parameters,
					int position_word, Vector<int>& ref) = 0;
    
    virtual void SetBoundaryConditionMesh(int, int) = 0;
    virtual void AddPeriodicConditionMesh(const TinyVector<int, 2>&, int) = 0;

    virtual void AddDirichletPointConditions(Vector<IVect>& ) = 0;
    
    
    /***********************
     * Dirichlet condition *
     ***********************/
    
  public :
    // Inline methods
    int GetNbDirichletDof() const;
    int GetNbGlobalDirichletDof() const;
    int GetDirichletDofNumber(int) const;
    bool IsDofDirichlet(int i) const;
    const Vector<bool>& GetIsDofDirichlet() const;
    bool UseSymmetricDirichlet() const;
    void EnableSymmetricDirichlet(bool sym = true);
    Real_wp GetCoefficientDirichlet() const;
    void SetCoefficientDirichlet(const Real_wp& coef);
    void SetDirichletDof(int i, bool b);
    const IVect& GetDirichletDofNumber() const;
    
    int GetNbSupportedComponents(int) const;
    int GetSupportedComponent(int, int) const;
    void SetSupportedComponents(int, const IVect&);
    
    template<class T>
    void ImposeNullDirichletCondition(Vector<T>& xsol) const;

    int GetHighConductivityOrder() const;

    // other methods related to Dirichlet
    template<class Dimension>
    void FindDofsOnReference(const VarProblem<Dimension>& var,
			     const Vector<int>& ref_cond, int ref_target, Vector<int>& Dofs);

  protected:
    template<class Dimension>
    void FindLocalDirichletDof(int i, int num_face, int ne, int num_loc, int ref_cond,
                               const VarProblem<Dimension>& var,
                               Vector<int>& IndexDirichlet,
			       const ElementReference_Dim<Dimension>& Fb, int num = 0);

  public:
    template<class Dimension>
    void TreatDirichletCondition(const VarProblem<Dimension>& var);
    
    void SetDirichletDofs(int, const IVect&);
    void ResizeNbDof(int n);
    
    template<class T>
    Real_wp ComputeDirichletCoef(VirtualMatrix<T>& A);

    void ApplyDirichletCondition(const SeldonTranspose&, FemMatrixFreeClass_Base<Real_wp>& mat_sp,
				 Vector<Real_wp>& b_rhs, int k = 0) const;

    void ApplyDirichletCondition(const SeldonTranspose&, FemMatrixFreeClass_Base<Real_wp>& mat_sp,
				 Vector<Complex_wp>& b_rhs, int k = 0) const;

    void ApplyDirichletCondition(const SeldonTranspose&, FemMatrixFreeClass_Base<Complex_wp>& mat_sp,
				 Vector<Real_wp>& b_rhs, int k = 0) const;

    void ApplyDirichletCondition(const SeldonTranspose&, FemMatrixFreeClass_Base<Complex_wp>& mat_sp,
				 Vector<Complex_wp>& b_rhs, int k = 0) const;
        
    void UpdateDirichletDofs();


    /**********************
     * Periodic condition *
     **********************/
    
    // Inline methods
    virtual int GetNbModes() const;
    virtual int GetNbModesSource() const;
    int GetModeNumber(int n) const;
    int GetCurrentModeNumber() const;
    void SetCurrentModeNumber(int n);
    bool ModesNotStored() const;
    void ForceStorageModes(bool store = true);
    int GetSymmetryType() const;
    int GetNbPeriodicDof(int n = 0) const;
    int GetPeriodicDof(int, int n = 0) const;
    int GetOriginalPeriodicDof(int, int n = 0) const;
    int GetFormulationForPeriodicCondition() const;
    
    void SetModesToCompute(const IVect& num);
    void PushBackMode(int n);
    
    int GetNbPeriodicModesX() const;
    int GetNbPeriodicModesY() const;
    int GetNbPeriodicModesZ() const;
    
    // other methods
    void GetPeriodicNumberModes(int& nx, int& ny, int& nz, bool& teta_sym) const;
    void GetPeriodicModes(int n, int& kx, int& ky, int& kz, bool& teta_sym) const;
    void GetPeriodicModes(int n, Complex_wp& kx, Complex_wp& ky, Complex_wp& kz) const;
    
    // modifies matrix/rhs due to periodic condition
    template<class T, class Storage>
    void SetPeriodicCondition(Matrix<T, Symmetric, Storage>&);

    template<class T>
    void SetPeriodicCondition(Matrix<T, Symmetric, DiagonalRow>&);
    
    template<class T>
    void SetPeriodicCondition(Matrix<T, General, BlockDiagRow>&);
    
    template<class T>
    void SetPeriodicCondition(FemMatrixFreeClass_Base<T>&);
    
    template<class T>
    void SetPeriodicCondition(Matrix<T, General, ArrayRowSparse>&);
    
    void SetPeriodicCondition(Matrix<Complex_wp, General, ArrayRowComplexSparse>&);
    
#ifdef SELDON_WITH_MPI
    template<class T>
    void SetPeriodicCondition(DistributedMatrix<T, General, ArrayRowSparse>&);

    void SetPeriodicCondition(DistributedMatrix<Complex_wp, General, ArrayRowComplexSparse>&);    
#endif

    template<class T>
    void ApplyPeriodicCondition(Vector<T>& brhs) const;
    

    /*******************************
     * Absorbing condition and PML *
     *******************************/
    
    // Inline methods
    int GetOrderAbsorbingCondition() const;
    int GetNbEltPML() const;
    int GetNbGlobalEltPML() const;
    //Real_wp GetOmegaPML() const;
    //void SetOmegaPML(const Real_wp&);
    
    virtual void MltMuIntegrationByParts(int ref, int ne, int num_loc, int k, Real_wp&) const;
    virtual void MltMuIntegrationByParts(int ref, int ne, int num_loc, int k, Complex_wp&) const;

    virtual void SetPhysicalIndexAtInfinity(const Vector<bool>&);
    const Real_wp& GetMaximumVelocityPML() const;    


    // Methods related to PML
    template<class Dimension>
    void FindElementsInsidePML(DistributedProblem<Dimension>& var);
    
    void EvaluateFunctionTauPML(const Real_wp& dx, const Real_wp& dsig, const Real_wp& a,
                                Real_wp& zeta, Real_wp& zeta_prime) const;

    template<class Dimension>
    void GetDampingFactorPML(Mesh<Dimension>& mesh, int num_pml, int type_pml,
			     const typename Dimension::R_N& point,
			     typename Dimension::R_N_Complex_wp& zeta,
			     typename Dimension::R_N_Complex_wp& point_prime);
    
    void GetDampingTauPML(const Mesh<Dimension2>& mesh, int num_pml, int type_pml,
			  const R2& point, R2& zeta, R2& zeta_prime) const;

    void GetDampingTauPML(const Mesh<Dimension3>& mesh, int num_pml, int type_pml,
			  const R3& point, R3& zeta, R3& zeta_prime) const;
    
    // other methods    
    template<class T, class Dimension>
    void AddMatrixImpedanceBoundary(const T& alpha, const Vector<int>& ref_cond, int ref_target,
				    const GlobalGenericMatrix<T>& nat_mat,
				    VirtualMatrix<T>& mat_sp, int offset_row, int offset_col,
				    ImpedanceFunction_Base<T, Dimension>& fimpedance,
				    bool changeColonne, bool changeLigne,
				    const VarProblem<Dimension>& var);
    
    void CopyInputData(const VarBoundaryCondition_Base& var);

    virtual Complex_wp ComputeImpedanceCoefficient(const IVect& ref, VectComplex_wp& U);
    
    // pure virtual methods
    virtual void AddBoundaryConditionTerms(VirtualMatrix<Real_wp>&,
					   const GlobalGenericMatrix<Real_wp>&,
					   int offset_row = 0, int offset_col = 0) = 0;
    
    virtual void AddBoundaryConditionTerms(VirtualMatrix<Complex_wp>&,
					   const GlobalGenericMatrix<Complex_wp>&,
					   int offset_row = 0, int offset_col = 0) = 0;

    virtual void InitCyclicDomain() = 0;
    virtual void TreatDirichletCondition() = 0;
    virtual void FindElementsInsidePML() = 0;
    virtual void ComputeQuasiPeriodicPhase() = 0;
    virtual void AllocateTauPML() = 0;

    virtual void GetPeriodicDofNumbers(int, int&, int&, int n = 0) const = 0;
    virtual void GetPeriodicPhase(int i, Real_wp&) const = 0;
    virtual void GetPeriodicPhase(int i, Complex_wp&) const = 0;
    
    virtual void MltParamCondition(int ref, int k, Real_wp&) const = 0;
    virtual void MltParamCondition(int ref, int k, Complex_wp&) const = 0;
    
  };


  //! Specificity of boundary conditions for real/complex numbers and/or 2-D/3-D
  template<class Complexe, class Dimension>
  class VarBoundaryCondition_Dim : public VarBoundaryCondition_Base
  {
  public :
    typedef typename Dimension::R_N R_N;
    typedef typename Dimension::R_N_Complex_wp R_N_Complex_wp;
    typedef Vector<R_N> VectR_N;
    typedef TinyVector<Complexe, Dimension::dim_N> R_N_Complexe;
    
  protected :
    //! if the user wants points with Dirichlet condition
    Vector<R_N> point_dirichlet_condition;

    //! parameters of the boundary condition (impedance by example)
    Vector<Vector<Complexe> > Param_condition;

    Vector<Complexe> periodic_dof_phase;

    Complexe Coef_Imped_ABC; //!< impedance coefficient for ABC

    //! tau damping for PML layers
    Vector<Vector<R_N_Complexe> > Glob_TauPML;
    Vector<Vector<R_N_Complexe> > Glob_TauPrimePML; //!< derivative of tau damping for PML layers
    
  private :
    //! reference to the mesh object
    Mesh<Dimension>& mesh_reference;
    DistributedProblem<Dimension>& var_comm;

    void InitDefaultValues();

  public :
    template<class TypeEquation>
    VarBoundaryCondition_Dim(EllipticProblem<TypeEquation>&);
    
    // Inline methods
    void GetPeriodicPhase(int i, Real_wp&) const;
    void GetPeriodicPhase(int i, Complex_wp&) const;

    void GetPeriodicDofNumbers(int i, int& num_per, int& num_orig, int n = 0) const;
    
    int GetProcOriginalPeriodicDof(int i) const;

    const Complexe& GetImpedanceCoefficientABC() const;

    const TinyVector<Complexe, Dimension::dim_N>& GetTauPML(int, int) const;
    const Complexe& GetTauPML(int, int, int) const;
    const Complexe& GetPrimitiveTauPML(int, int, int) const;

    void GetDampingFactorPML(const R_N& point, int num_pml, 
                             int type_pml, R_N& zeta, R_N& zeta_prime);

    Vector<Vector<Complexe> >& GetParamCondition();
    const Complexe& GetParamCondition(int, int) const;
    virtual void MltParamCondition(int ref, int k, Real_wp&) const;
    virtual void MltParamCondition(int ref, int k, Complex_wp&) const;

    void SetBoundaryConditionMesh(int ref, int type);
    void AddPeriodicConditionMesh(const TinyVector<int, 2>& ref, int type);

    void GetDampingFactorPML(const R_N& point, int num_pml, int type_pml,
                             R_N_Complex_wp& zeta, R_N_Complex_wp& zeta_prime);
    
    // Other methods
    void GetMemoryUsed(map<string, size_t>& var) const;

    const Matrix<int>& GetGeometryPhaseData(int iquad, int num_loc,
					    int& nf, int& ref_boundary, int& ref, bool& new_face,
					    int& rot, Complexe& phase, Complexe& phase_conj) const;
    
    void SetInputData(const string& description_field, const VectString& parameters);

  protected :
    void CheckAndReallocateReferences(int ref);
    void FillParameterCondition(const Vector<string>& parameters,
				int position_word, Vector<int>& ref);
    
    void AddDirichletPointConditions(Vector<IVect>& IndexOther);
    
  public :
    void AllocateTauPML();
    void ComputeStoreCoefficientsPML(int i1, int i, const VectR_N& AllPoints);
    
    void InitCyclicDomain();    

    void ComputeQuasiPeriodicPhase();

    void CopyInputData(const VarBoundaryCondition_Base& var);
    
  };

  
  template<class Complexe, class Dimension>
  class VarBoundaryCondition : public VarBoundaryCondition_Dim<Complexe, Dimension>
  {
  private:
    VarProblem<Dimension>& var_problem;
#ifdef MONTJOIE_WITH_TRANSMISSION
    VarTransmission_Base<Dimension>& var_transmission_base;
#endif
    VarGeneralizedImpedance_Base& var_gibc_base;
    
  public:
    template<class TypeEquation>
    VarBoundaryCondition(EllipticProblem<TypeEquation>& var);
    
    // Inline methods
    void FindDofsOnReference(const Vector<int>& ref, Vector<int>& Dofs);
    void TreatDirichletCondition();

#ifdef MONTJOIE_WITH_TRANSMISSION
    const VarTransmission_Base<Dimension>& GetTransmissionProblem() const;
#endif
    
    template<class T>
    void AddMatrixImpedanceBoundary(const T& alpha, const Vector<int>& ref_cond, int ref_target,
				    const GlobalGenericMatrix<T>& nat_mat,
				    VirtualMatrix<T>& mat_sp, int offset_row, int offset_col,
                                    ImpedanceFunction_Base<T, Dimension>& fimpedance,
				    bool changeColonne = false, bool changeLigne = false);
    
    void FindElementsInsidePML();
        
    void AddBoundaryConditionTerms(VirtualMatrix<Real_wp>&,
				   const GlobalGenericMatrix<Real_wp>&,
				   int offset_row = 0, int offset_col = 0);

    void AddBoundaryConditionTerms(VirtualMatrix<Complex_wp>&,
				   const GlobalGenericMatrix<Complex_wp>&,
				   int offset_row = 0, int offset_col = 0);
    
    // General method to add boundary terms in finite element matrix
    template<class T>
    void AddBoundaryConditionGen(VirtualMatrix<T>& mat_sp,
				 const GlobalGenericMatrix<T>& nat_mat,
				 int offset_row = 0, int offset_col = 0);

    // pure virtual methods to allocate impedance objects
    virtual ImpedanceFunction_Base<Real_wp, Dimension>*
    GetNewImpedanceABC(const Real_wp&) const = 0;
    
    virtual ImpedanceFunction_Base<Complex_wp, Dimension>*
    GetNewImpedanceABC(const Complex_wp&) const = 0;

    virtual ImpedanceFunction_Base<Real_wp, Dimension>*
    GetNewImpedanceGeneric(const Real_wp&) const = 0;
    
    virtual ImpedanceFunction_Base<Complex_wp, Dimension>*
    GetNewImpedanceGeneric(const Complex_wp&) const = 0;
    
#ifdef MONTJOIE_WITH_HIGH_CONDUCTIVITY_MODEL    
    virtual ImpedanceFunction_Base<Real_wp, Dimension>*
    GetNewImpedanceHighConductivity(const Real_wp&) const = 0;
    
    virtual ImpedanceFunction_Base<Complex_wp, Dimension>*
    GetNewImpedanceHighConductivity(const Complex_wp&) const = 0;
#endif

  };

  
  //! General class for impedance
  template<class T, class Dimension>
  class VirtualImpedanceFunction
  {
  protected:
    Vector<Vector<T> > stored_coef_phi; //!< coefficient for a scalar function f
    Vector<Vector<T> > stored_coef_grad; //!< coefficient for a scalar function g
    Vector<T> coef_phi, coef_grad;
    bool read_param_condition;
    bool only_surface_dof_involved;
    Vector<typename Dimension::R_N> vec_normale; //!< normale for each quadrature point
    
  private:    
    const DistributedProblem<Dimension>& var_problem;
    const VarBoundaryCondition_Base& var_boundary;
    
    void InitDefaultValues();
    
  public:
    VirtualImpedanceFunction(const DistributedProblem<Dimension>&);
    
    virtual ~VirtualImpedanceFunction();
    
    // Inline methods
    bool InvolveOnlyTangentialDofs() const;
    void SetCoefficient(const T& a, const T& b);
    
    inline T GetStoredCoefficient(int i, int j) const { return stored_coef_phi(i)(j); }
    inline T GetStoredGradCoefficient(int i, int j) const { return stored_coef_grad(i)(j); }
    
    virtual bool PresenceGradient() const;
    virtual T GetCoefficient(int i, int iquad, int num_loc, int k, int ref_d, int ref,
			     const SetPoints<Dimension>&, const SetMatrices<Dimension>&);
    
    // method that will compute and store f
    virtual void EvaluateImpedancePhi(int i, int iquad, int num_edge, int num_loc, int k,
				      const GlobalGenericMatrix<T>& nat_mat, int ref_domain,
				      const SetPoints<Dimension>&, const SetMatrices<Dimension>&);
    
    // method that will compute and store g
    virtual void EvaluateImpedanceGrad(int i, int iquad, int num_edge, int num_loc, int k,
				       const GlobalGenericMatrix<T>& nat_mat, int ref_domain,
				       const SetPoints<Dimension>&, const SetMatrices<Dimension>&);
    
  };
  
  
  //! impedance function (diagonal coefficients)
  template<class T>
  class ImpedanceFunction_Base<T, Dimension2>
    : public VirtualImpedanceFunction<T, Dimension2>
  {
  private:    
    const VarProblem<Dimension2>& var_problem;
    const VarBoundaryCondition_Base& var_boundary;
    
    void InitDefaultValues();
    
  public :    
    ImpedanceFunction_Base(const VarProblem<Dimension2>&);    
    
    // computes f(phi_i, dphi_i) (m = unknown number associated with row i) for H^1 unknowns
    virtual void ApplyImpedancePhi_H1(int m, int j, int, const TinyVector<Real_wp, 1>& phi,
                                      const R2& grad_phi, Vector<T>& f_phi);

    // computes g(phi_i, dphi_i) (m = unknown number associated with row i) for H^1 unknowns
    virtual void ApplyImpedanceGrad(int m, int j, int, const TinyVector<Real_wp, 1>& phi,
                                       const R2& grad_phi, Vector<T>& g_phi);
    
    // computes f(phi_i, dphi_i) (m = unknown number associated with row i) for H(curl) unknowns
    virtual void ApplyImpedancePhi_Hcurl(int m, int j, int, const TinyVector<Real_wp, 2>& phi,
                                         const TinyVector<Real_wp, 1>& curl_phi, Vector<T>& f_phi);

    // computes g(phi_i, dphi_i) (m = unknown number associated with row i) for H(curl) unknowns
    virtual void ApplyImpedanceCurl(int m, int j, int, const TinyVector<Real_wp, 2>& phi,
				    const TinyVector<Real_wp, 1>& curl_phi, Vector<T>& g_phi);

    // computes f(phi_i, dphi_i) (m = unknown number associated with row i) for H(div)
    virtual void ApplyImpedancePhi_Hdiv(int m, int j, int, const R2& phi,
				   const TinyVector<Real_wp, 1>& div_phi, Vector<T>& f_phi);
    
    // computes g(phi_i, dphi_i) (m = unknown number associated with row i) for H(div)
    virtual void ApplyImpedanceDiv(int m, int j, int, const R2& phi,
                                   const TinyVector<Real_wp, 1>& div_phi, Vector<T>& g_phi);
  };


  //! impedance function (diagonal coefficients)
  template<class T>
  class ImpedanceFunction_Base<T, Dimension3>
    : public VirtualImpedanceFunction<T, Dimension3>
  {
  private:    
    const VarProblem<Dimension3>& var_problem;
    const VarBoundaryCondition_Base& var_boundary;
    
    void InitDefaultValues();
    
  public :    
    ImpedanceFunction_Base(const VarProblem<Dimension3>&);    
    
    // computes f(phi_i, dphi_i) (m = unknown number associated with row i) for H^1 unknowns
    virtual void ApplyImpedancePhi_H1(int m, int j, int, const TinyVector<Real_wp, 1>& phi,
                                      const R3& grad_phi, Vector<T>& f_phi);

    // computes g(phi_i, dphi_i) (m = unknown number associated with row i) for H^1 unknowns
    virtual void ApplyImpedanceGrad(int m, int j, int, const TinyVector<Real_wp, 1>& phi,
				   const R3& grad_phi, Vector<T>& g_phi);
    
    // computes f(phi_i, dphi_i) (m = unknown number associated with row i) for H(curl)
    virtual void ApplyImpedancePhi_Hcurl(int m, int j, int, const TinyVector<Real_wp, 3>& phi,
				   const TinyVector<Real_wp, 3>& curl_phi, Vector<T>& f_phi);

    // computes g(phi_i, dphi_i) (m = unknown number associated with row i) for H(curl)
    virtual void ApplyImpedanceCurl(int m, int j, int, const TinyVector<Real_wp, 3>& phi,
				    const TinyVector<Real_wp, 3>& curl_phi, Vector<T>& g_phi);
    
    // computes f(phi_i, dphi_i) (m = unknown number associated with row i) for H(div)
    virtual void ApplyImpedancePhi_Hdiv(int m, int j, int, const R3& phi,
                                        const TinyVector<Real_wp, 1>& div_phi, Vector<T>& f_phi);
    
    // computes g(phi_i, dphi_i) (m = unknown number associated with row i) for H(div)
    virtual void ApplyImpedanceDiv(int m, int j, int, const R3& phi,
                                   const TinyVector<Real_wp, 1>& div_phi, Vector<T>& g_phi);
    
  };

  
  //! class for impedance boundary condition (ConditionReference = ref IMPEDANCE)
  template<class T, class TypeEquation>
  class ImpedanceGeneric
    : public ImpedanceFunction_Base<T, typename TypeEquation::Dimension>
  {
  public:
    ImpedanceGeneric(const EllipticProblem<TypeEquation>& var);

  };
  
  
  //! class for absorbing boundary condition (ConditionReference = ref ABSORBING)
  template<class T, class TypeEquation>
  class ImpedanceABC
    : public ImpedanceFunction_Base<T, typename TypeEquation::Dimension>
  {
  public :
    ImpedanceABC(const EllipticProblem<TypeEquation>& var);

  };
  
  
#ifdef MONTJOIE_WITH_HIGH_CONDUCTIVITY_MODEL
  //! empty class, overloaded
  template<class T, class TypeEquation>
  class ImpedanceHighConductivity
    : public ImpedanceFunction_Base<T, typename TypeEquation::Dimension>
  {
  public:
    ImpedanceHighConductivity(const EllipticProblem<TypeEquation>& var);
    
  };
#endif
  
} // namespace Montjoie

#define MONTJOIE_FILE_BOUNDARY_CONDITION_HARMONIC_HXX
#endif

