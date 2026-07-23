#ifndef MONTJOIE_FILE_MAXWELL_AXISYM_HARMONIC_HXX

// headers for resolution of axi-symmetric maxwell equations
namespace Montjoie
{
  
  //! class to solve time-harmonic Maxwell equations on axisymmetric domains
  template<class Complexe>
  class HarmonicMaxwell_Axi : public VarAxisymProblem
  {
  public :
    // typedef declarations
    typedef All_MatrixLU<Complexe> MatrixLU; //!< Lu factors
    typedef Vector<Complexe> VectComplexe; //!< vector of floats or complex floats

    //! dielectric permittivity        
    Vector<TensorPhysicalIndice<Dimension2, 3, Complexe> > ref_epsilon; 
    //! inverse of dielectric permittivity
    Vector<TensorPhysicalIndice<Dimension2, 3, Complexe> > ref_invEpsilon; 
    //! magnetic permeability
    Vector<TensorPhysicalIndice<Dimension2, 3, Complexe> > ref_mu; 
    //! inverse of magnetic permeability
    Vector<TensorPhysicalIndice<Dimension2, 3, Complexe> > ref_invMu; 
    //! damping (scalar)
    Vector<ScalarPhysicalIndice<Dimension2, Complexe> > ref_sigma;
    
    //! radius on each quadrature point (\tilde{r} in PML domain)
    Vector<VectComplexe> Glob_rtilde;
 
    //! evaluation of J on a surface
    string name_file_outputJ;
    int nb_subdiv_outputJ; bool output_Hy_teta;
    IVect ref_outputJ;

    //! dieletric permittivity in vacuum
    Complexe epsilon0;
    //! magnetic permeability in vacuum
    Complexe mu0, invMu0;
    
  private:
    DistributedProblem<Dimension2>& var_problem;
    VarBoundaryCondition_Dim<Complexe, Dimension2>& var_boundary;

  public:
    //! default constructor
    template<class TypeEquation>
    HarmonicMaxwell_Axi(EllipticProblem<TypeEquation>& var);
        
    // Inline methods
    int GetNbPhysicalIndices() const;
    bool IsVaryingMedia(int i) const;

    Real_wp GetVelocityOfMedia(int ref) const;
    Real_wp GetVelocityOfInfinity() const;

    void GetVaryingIndices(Vector<PhysicalVaryingMedia<Dimension2, Complex_wp>* >& rho_complex,
			   Vector<PhysicalVaryingMedia<Dimension2, Real_wp>* >& rho_real, IVect& num_ref,
			   IVect& num_index, IVect& num_component, Vector<bool>& compute_grad,
			   Vector<bool>& compute_hess);
        
    // additionnal input datas
    void SetInputData(const string& description_field, const VectString& parameters);
    
    // allocation of arrays containing physical indices
    void InitIndices(int n);
    
    // set physical indices of media with reference i
    void SetIndices(int i, const VectString& parameters);
    void SetPhysicalIndex(const string& name_media, int i, const VectString& parameters);
    string GetPhysicalIndexName(int m) const;
    
    void GetVaryingIndex(Vector<PhysicalVaryingMedia<Dimension2, Complexe>* >& rho_complex, 
			 IVect& num_ref, IVect& num_index, IVect& num_component,
			 Vector<bool>& compute_grad, Vector<bool>& compute_hess);
    
    void FinalizeComputationVaryingIndices();    

    void ModifyPMLCoefficient(TinyMatrix<Complexe, Symmetric, 2, 2>& epsilon_rz,
			      Complexe& epsilon_teta,
			      TinyMatrix<Complexe, Symmetric, 2, 2>& mu_rz,
			      Complexe& mu_teta, Complexe&, Complexe&, int i1, int j) const;
    
    void AllocateMassMatrices();
    void ComputeLocalMassMatrix(int i);
    
  };


  template<class TypeEquation>
  class HarmonicMaxwellAxi_Eq : public HarmonicMaxwell_Axi<typename TypeEquation::Complexe>,
				public VarHarmonic<TypeEquation>
  {
    typedef typename TypeEquation::Complexe Complexe;
    
  public:
    HarmonicMaxwellAxi_Eq();

    void SetWaveVector(const R3& k);
    const R3& GetWaveVector() const;
    const R3& GetPhaseOrigin() const;
    int GetNbModesSource() const;
    
    bool UseNumericalIntegration(int iquad) const;
    bool IsSymmetricProblem(bool eigen = false) const;
    
    int GetNbPhysicalIndices() const;
    bool IsVaryingMedia(int i) const;
    bool IsVaryingMedia(int m, int i) const;

    Real_wp GetVelocityOfMedia(int ref) const;
    Real_wp GetVelocityOfInfinity() const;
    
    void SetInputData(const string& description_field, const VectString& parameters);
    void CheckInputMesh();
    
    void InitIndices(int n);
    
    void SetIndices(int i, const VectString& parameters);
    void SetPhysicalIndex(const string& name_media, int i, const VectString& parameters);
    string GetPhysicalIndexName(int m) const;
    
    void FinalizeComputationVaryingIndices();    

    void AllocateMassMatrices();
    void ComputeLocalMassMatrix(int i, int N, bool linear_sparse,
				SetPoints<Dimension2>& PointsElem,
				SetMatrices<Dimension2>& MatricesElem,
				IVect& OrderFace, const ElementGeomReference<Dimension2>& Fb);
    
    void GetVaryingIndices(Vector<PhysicalVaryingMedia<Dimension2, Complex_wp>* >& rho_complex,
			   Vector<PhysicalVaryingMedia<Dimension2, Real_wp>* >& rho_real, IVect& num_ref,
			   IVect& num_index, IVect& num_component, Vector<bool>& compute_grad,
			   Vector<bool>& compute_hess);
    
  };
  
  
  //! class to describe time-harmonic Maxwell equations in axisymmetric equations
  class HarmonicMaxwellEquation_HcurlAxi : public GenericEquation<Complex_wp>
  {
  public :
    typedef Dimension2 Dimension;
    
    enum {nb_unknowns = 2, nb_unknowns_scal = 2, nb_unknowns_hdg=0,
	  nb_components_en = 2, nb_components_hn = 2, nb_unknowns_vec = 3,
          type_element = 2};
    
    static bool SymmetricGlobalMatrix();
    static bool SymmetricElementaryMatrix();

    static inline IVect GetOtherElementType() { IVect num(1); num(0) = 1; return num; }
    
    template<class TypeEquation>
    static void ComputeMassMatrix(EllipticProblem<TypeEquation>& var,
                                  int i, const ElementReference_Dim<Dimension2>& Fb);

    static Real_wp GetCoefficientMassMatrix(const VarProblem_Base& vars);
    static Real_wp GetCoefficientStiffnessMatrix(const VarProblem_Base& vars);
    static Real_wp GetCoefMode(const VarProblem_Base& vars);
    
  };


  //! class for Silver-Muller condition for Maxwell equations in axisymmetric domain
  template<>
  class ImpedanceABC<Complex_wp, HarmonicMaxwellEquation_HcurlAxi>
    : public ImpedanceFunction_Base<Complex_wp, Dimension2>
  {
  public :
    VectR2 VecNormale; VectReal_wp VecRadius;
    Vector<Complex_wp> VecImpedance;
    bool modified_formulation;
    int number_mode;

  private:
    const EllipticProblem<HarmonicMaxwellEquation_HcurlAxi>& var_maxwell;
    
  public:
    template<class TypeEquation>
    ImpedanceABC(const EllipticProblem<TypeEquation>&);

    bool PresenceGradient() const;
    
    void EvaluateImpedancePhi(int i, int num_elem, int num_edge, int num_loc, int k,
			      const GlobalGenericMatrix<Complex_wp>& nat_mat, int ref,
			      const SetPoints<Dimension2>& Pts,
			      const SetMatrices<Dimension2>& Mat);
    
    void ApplyImpedancePhi_H1(int n, int j, int offset,
                              const TinyVector<Real_wp, 1>& phi_boundary,
                              const TinyVector<Real_wp, 2>& dphi_boundary, VectComplex_wp& feval_phi);

    void ApplyImpedancePhi_Hcurl(int n, int j, int offset,
                                 const TinyVector<Real_wp, 2>& phi_boundary,
                                 const TinyVector<Real_wp, 1>& dphi_boundary, VectComplex_wp& feval_phi);

    void ApplyImpedanceGrad(int n, int j, int, const TinyVector<Real_wp, 1>& phi_boundary,
			    const TinyVector<Real_wp, 2>& dphi_boundary, VectComplex_wp& feval_dphi);

    void ApplyImpedanceCurl(int n, int j, int, const TinyVector<Real_wp, 2>& phi_boundary,
			    const TinyVector<Real_wp, 1>& dphi_boundary, VectComplex_wp& feval_dphi);

  };


  //! base class for impedance boundary condition
  template<>
  class ImpedanceGeneric<Complex_wp, HarmonicMaxwellEquation_HcurlAxi>
    : public ImpedanceFunction_Base<Complex_wp, Dimension2>
  {
  public :
    template<class TypeEquation>
    ImpedanceGeneric(const EllipticProblem<TypeEquation>&);

    void EvaluateImpedancePhi(int i, int num_elem, int num_edge, int num_loc, int k,
			      const GlobalGenericMatrix<Complex_wp>& nat_mat, int ref_d,
			      const SetPoints<Dimension2>& Pts, const SetMatrices<Dimension2>& Mat);

    void ApplyImpedancePhi_H1(int n, int j, int offset,
                              const TinyVector<Real_wp, 1>& phi_boundary,
                              const TinyVector<Real_wp, 2>& dphi_boundary, VectComplex_wp& feval_phi);

    void ApplyImpedancePhi_Hcurl(int n, int j, int offset,
                                 const TinyVector<Real_wp, 2>& phi_boundary,
                                 const TinyVector<Real_wp, 1>& dphi_boundary, VectComplex_wp& feval_phi);
    
  };
  
  
  //! class for the computation of radar cross section for axisymmetric Maxwell's equations
  class VarComputationRCS_MaxwellAxi : public VarComputationRCS_Axi
  {
  private:
    EllipticProblem<HarmonicMaxwellEquation_HcurlAxi>& var_maxwell;
    
  public :    
    template<class TypeEquation>
    VarComputationRCS_MaxwellAxi(EllipticProblem<TypeEquation>& var);

    void ComputeRCS(const VectReal_wp& U0);
    void ComputeRCS(const VectComplex_wp& U0);

    void ComputeRCS(const Matrix<Complex_wp, General, ColMajor>& U0, VectR3_Complex_wp& RCS_Vector);
        
    void ComputeEnHn(const VectComplex_wp& U0, VectR3_Complex_wp& traceEn,
		     VectR3_Complex_wp& trace_Hn);
    
    void AddContributionMode(int num_angle, const R3& ur, const VectR3_Complex_wp& trace_En,
			     const VectR3_Complex_wp& trace_Hn, R3_Complex_wp& rcs_teta);
    
    void ComputeMonostaticRCS(const Matrix<Complex_wp, General, ColMajor>& U0, VectR3_Complex_wp& rcs_vector);
    void WriteMonostaticRCS(const VectR3_Complex_wp& RCS_Vector);
    
  };
  
  
  //! class for computing radar cross section for axisymmetric Maxwell's equations
  template<>
  class VarComputationRCS<HarmonicMaxwellEquation_HcurlAxi>
    : public VarComputationRCS_MaxwellAxi
  {
  public :
    template<class TypeEquation>
    VarComputationRCS(EllipticProblem<TypeEquation>& var);
    
  };

  
  //! time-harmonic maxwell's equations in axisymmetric domain with LDG
  class HarmonicMaxwellEquationAxiDG : public GenericEquation<Complex_wp>
  {
  public :
    typedef Dimension2 Dimension;
    
    enum {nb_unknowns = 6, nb_components_en = 2, nb_components_hn = 2,nb_unknowns_hdg=0};
    
    static const bool FirstOrderFormulation = true;

    static inline bool SymmetricGlobalMatrix() { return true; }
    static inline bool SymmetricElementaryMatrix() { return false; }
    
    template<class TypeEquation, class T0, class MatMass>
    static void GetTensorMass(const EllipticProblem<TypeEquation>& vars,
			      int i, int j, const GlobalGenericMatrix<T0>& nat_mat, int ref, MatMass& mass);
    
    template<class TypeEquation, class T0, class Vector1>
    static void ApplyTensorMass(const EllipticProblem<TypeEquation>& var, int i, int j,
				const GlobalGenericMatrix<T0>& nat_mat, int ref, Vector1& U, Vector1& V);
    
    template<class TypeEquation, class T0, class Vector1>
    static void GetNeededDerivative(const EllipticProblem<TypeEquation>& vars,
                                    const GlobalGenericMatrix<T0>& nat_mat,
                                    Vector1& unknown_to_derive, Vector1& fct_test_to_derive);
    
    template<class TypeEquation, class T0, class MatStiff>
    static void GetGradPhiTensor(const EllipticProblem<TypeEquation>& vars,
				 int num_elem, int jloc, const GlobalGenericMatrix<T0>& nat_mat, int ref,
				 MatStiff& Dgrad_phi, MatStiff& Ephi_grad);
    
    template<class TypeEquation, class T0, class Vector1, class Vector2>
    static void ApplyGradientUnknown(const EllipticProblem<TypeEquation>& var,
				     int i, int j, const GlobalGenericMatrix<T0>& nat_mat,
                                     int ref, Vector1& dU, Vector2& V);

    template<class TypeEquation, class T0, class Vector1, class Vector2>
    static void ApplyGradientFctTest(const EllipticProblem<TypeEquation>& var,
				     int i, int j, const GlobalGenericMatrix<T0>& nat_mat,
                                     int ref, Vector1& U, Vector2& dV);
    
    template<class Matrix1, class GenericPb, class T0>
    static void GetNabc(Matrix1& Nabc, R2& normale,
			int ref, int iquad, int k, const GlobalGenericMatrix<T0>& nat_mat, int ref_d,
			const GenericPb& vars, const ElementReference<Dimension2, 1>& Fb);

    template<class Vector1, class TypeEquation, class T0>
    static void MltNabc(R2& normale, int ref,
                        const Vector1& Vn, Vector1& Un, int num_elem1, int k,
                        const GlobalGenericMatrix<T0>& nat_mat, int ref_d,
                        const EllipticProblem<TypeEquation>& vars, const ElementReference<Dimension2, 1>& Fb);
    
    template<class Matrix1, class GenericPb, class T0>
    static void GetPenalDG(Matrix1& Nabc, R2& normale,
			   int iquad, int k, int nf, const GlobalGenericMatrix<T0>& nat_mat, int ref, int ref2,
			   const GenericPb& vars, const ElementReference<Dimension2, 1>& Fb);
    
    template<class Vector1, class Vector2, class GenericPb, class T0>
    static void MltPenalDG(const R2& normale, const Vector1& Vn, Vector2& Un,
			   int i, int k, int nf, const GlobalGenericMatrix<T0>& nat_mat, int ref, int ref2,
			   const GenericPb& vars, const ElementReference<Dimension2, 1>& Fb);
    
  };
  
  
  //! class for the solution of axisymmetric Maxwell's equations with LDG formulation
  template<>
  class EllipticProblem<HarmonicMaxwellEquationAxiDG> 
    : public HarmonicMaxwellAxi_Eq<HarmonicMaxwellEquationAxiDG>
  {
  public:
    
  };
  
}

#define MONTJOIE_FILE_MAXWELL_AXISYM_HARMONIC_HXX
#endif
