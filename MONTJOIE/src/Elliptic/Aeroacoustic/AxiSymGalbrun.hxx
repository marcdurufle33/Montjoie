#ifndef MONTJOIE_FILE_AXISYM_GALBRUN_HXX

namespace Montjoie
{

  class VarGalbrunIndex_Axi : public VarAxisymProblem
  {
  public :
    Vector<ScalarPhysicalIndice<Dimension2, Real_wp> > ref_rho0, ref_c0, ref_p0, ref_sigma;
    Vector<VectorPhysicalIndice<Dimension2, 3, Real_wp> > ref_m0, ref_g0; // AJOUT NATHAN
    
    Vector<VectR3> eval_flow;
    Vector<VectR2> grad_p0, grad_sigma, grad_rho, grad_gamma;
    Vector<Vector<TinyMatrix<Real_wp, General, 3, 2> > > grad_flow;
    Vector<Vector<TinyMatrix<Real_wp, General, 3, 2> > > grad_gravity;// AJOUT NATHAN
    Vector<VectReal_wp> eval_c0, eval_sigma, eval_rho, div_flow, eval_gamma;
    Vector<VectR3> eval_gravity; //AJOUT NATHAN

    enum {LEE_MODEL, SIMPLIFIED_LEE, BOGEY_BAILLY_JUVE, CONSERVATIVE, GALBRUN, GALBRUN_DIV};
    int type_model;
    bool apply_convective_derivate_source;
    bool compute_gravity; // AJOUT NATHAN

    int drop_unstable_terms;
    Real_wp coef_convective_term;
    enum { DROP_NONE, DROP_CONVECTIVE, DROP_NON_UNIFORM};

  private:
    VarProblem<Dimension2>& var_problem;
    VarOutputProblem_Dim<Dimension2>& var_output;
    
  public:
    template<class TypeEquation>
    VarGalbrunIndex_Axi(EllipticProblem<TypeEquation>&);
    
    // additionnal input datas
    void SetInputData(const string& description_field, const VectString& parameters);

    // allocation of arrays containing physical indices
    void InitIndices(int n);

    void GetVaryingIndices(Vector<PhysicalVaryingMedia<Dimension2, Real_wp>* >& rho_real, IVect& num_ref,
			   IVect& num_index, IVect& num_component, Vector<bool>& compute_grad,
			   Vector<bool>& compute_hess);

    void ComputePhysicalCoefficients();
    
    // set physical indices of media with reference i
    void SetIndices(int i, const VectString& parameters);
    void SetPhysicalIndex(const string& name_media, int i, const VectString& parameters);
    string GetPhysicalIndexName(int m) const;
    int GetNbPhysicalIndices() const;
    bool IsVaryingMedia(int i) const;
    Real_wp GetVelocityOfMedia(int ref) const;
    Real_wp GetVelocityOfInfinity() const;
    

    void ConstructAll(const string& name_file, const string& name_element,
		      Vector<string>& lines_data_file);
    
  };
  
  
  //! class to solve Galbrun's equation on axisymmetric domains
  template<class Complexe>
  class VarGalbrun_Axi : public VarGalbrunIndex_Axi
  {
  public :
    Vector<Vector<Complexe> > Glob_rtilde;
    Vector<Vector<Real_wp> > Glob_radius;

  private:
    DistributedProblem<Dimension2>& var_problem;
    VarBoundaryCondition_Dim<Complexe, Dimension2>& var_boundary;
    
  public:
    template<class TypeEquation>
    VarGalbrun_Axi(EllipticProblem<TypeEquation>&);
    
    void AllocateMassMatrices();    
    void ComputeLocalMassMatrix(int i);
    
  };
  
  
  //! volumetric source
  class VolumetricSource_AxiGalbrun
    : public VirtualSourceFEM<Complex_wp, Dimension2>
  {
  private:
    const DistributedProblem<Dimension2>& var_problem;
    const VarGalbrun_Axi<Complex_wp>& var_galbrun;
    const VarBoundaryCondition_Dim<Complex_wp, Dimension2>& var_boundary;

    GaussianSource<Dimension3> fsrc;    
    Vector<Complex_wp> polar;
    
  public :
    Complex_wp coef_vol;
    
    //! constructor with given problem
    template<class TypeEquation>
    VolumetricSource_AxiGalbrun(const EllipticProblem<TypeEquation>& var,
				const Vector<VectString>& param); 
    
    //! volumetric source
    bool IsNonNullVolumetricSource(const VectR2& s);
       
    //! evaluation of volumetric source
    void EvaluateVolumetricSource(int i, int j, const R2& x, Vector<Complex_wp>& f);
    
  };
  
  
  //! Galbrun's equation solved with LDG in cylindrical coordinates
  /*!
    rho_0 (-i omega + sigma + m_r d/dr - i m m_theta/r + m_z d/dz) u - rho_0 v = 0
    rho_0 (-i omega + sigma + m_r d/dr - i m m_theta/r + m_z d/dz) v_r - d/dr(rho_0 c_0^2 p)
         + (u_r/r -i m u_theta/r + du_z/dz) dp0/dr - du_z/dr dp0/dz = 0

    rho_0 (-i omega + sigma + m_r d/dr - i m m_theta/r + m_z d/dz) v_theta + i m/r rho_0 c_0^2 p
        + (i m u_r/r dp0/dr + i m u_z/r dp0/dz) = 0

    rho_0 (-i omega + sigma + m_r d/dr - i m m_theta/r + m_z d/dz) v_z - d/dz(rho_0 c_0^2 p)
        + (u_r/r -i m u_theta/r + du_r/dr) dp0/dz - du_r/dz dp0/dr = 0
    
    p - 1/r d/dr(r u_r) + i m/r u_theta - du_z/dz = 0
   */
  template<class T>
  class GalbrunAxiEquation_Base : public GenericEquation<T>
  {
  public :
    typedef Dimension2 Dimension;
    typedef GalbrunAxiEquation_Base<Real_wp> TypeEquationTime;
    
    static const bool FirstOrderFormulation = true;
    
    enum {nb_unknowns = 7, nb_components_u = 1, nb_unknowns_hdg=0,
          nb_components_grad = 2,
	  nb_components_en = 1, nb_components_hn = 1,
	  nb_unknowns_scal = 7, nb_unknowns_vec = 0};

    static inline bool SymmetricGlobalMatrix() { return false; }
    static inline bool SymmetricElementaryMatrix() { return false; }
    
    // for compatbility purpose
    template<class TypeEquation>
    static void ComputeMassMatrix(EllipticProblem<TypeEquation>& var,
                                  int num_elem, const ElementReference_Dim<Dimension>& Fb);
    
    // For a detailed description of the following methods
    // look at the file GenericEquation.hxx (class GenericEquation_Base)
    template<class GenericPb, class T0, class Vector1>
    static void GetNeededDerivative(const GenericPb& vars,
                                    const GlobalGenericMatrix<T0>& nat_mat,
				    Vector1& unknown_to_derive, Vector1& fct_test_to_derive);
    
    template<class TypeEquation, class T0, class Vector1>
    static void ApplyTensorStiffness(const EllipticProblem<TypeEquation>& var,
                                     int i, int j, const GlobalGenericMatrix<T0>& nat_mat,
                                     int ref, Vector1& dU, Vector1& dV);
    
    template<class TypeEquation, class T0, class MatStiff>
    static void GetGradGradTensor(const EllipticProblem<TypeEquation>& vars,
                                  int num_elem, int jloc, const GlobalGenericMatrix<T0>& nat_mat,
                                  int ref, MatStiff& Cgrad_grad);
    
    template<class TypeEquation, class T0, class MatStiff>
    static void GetGradPhiTensor(const EllipticProblem<TypeEquation>& vars,
                                 int num_elem, int jloc, const GlobalGenericMatrix<T0>& nat_mat, int ref,
                                 MatStiff& Dgrad_phi, MatStiff& Ephi_grad);
    
    template<class TypeEquation, class T0, class Vector1, class Vector2>
    static void ApplyGradientUnknown(const EllipticProblem<TypeEquation>& var,
				     int i, int j, const GlobalGenericMatrix<T0>& nat_mat,
                                     int ref, Vector1& dU, Vector2& V);

    template<class TypeEquation, class T0, class MatMass>
    static void GetTensorMass(const EllipticProblem<TypeEquation>& vars,
			      int i, int j, const GlobalGenericMatrix<T0>& nat_mat, int ref, MatMass& mass);
    
    template<class TypeEquation, class T0, class Vector1>
    static void ApplyTensorMass(const EllipticProblem<TypeEquation>& var, int i, int j,
				const GlobalGenericMatrix<T0>& nat_mat, int ref, Vector1& U, Vector1& V);

    template<class T0>
    static void GetAbsoluteD(TinyMatrix<T0, General, 7, 7>& Dtest, const R2& normale,
                             const Real_wp& c0, const T0& rtilde, bool non_unif, const Real_wp& alpha1_,
                             const Real_wp& alpha, const Real_wp& gamma, const Real_wp& br, const Real_wp& bz);
    
    template<class Matrix1, class GenericPb, class T0>
    static void GetAbsoluteD(Matrix1& Nabc, const R2& normale, int iquad, int k,
                             const GlobalGenericMatrix<T0>& nat_mat, const GenericPb& vars);
    
    template<class Matrix1, class GenericPb, class T0>
    static void GetNabc(Matrix1& Nabc, const R2& normale,
			int ref, int iquad, int k, const GlobalGenericMatrix<T0>& nat_mat, int ref_d, 
			const GenericPb& vars, const ElementReference<Dimension, 1>& Fb);
    
    template<class Vector1, class TypeEquation, class T0>
    static void MltNabc(const R2& normale, int ref,
                        const Vector1& Vn, Vector1& Un, int num_elem1,
                        int num_point, const GlobalGenericMatrix<T0>& nat_mat, int ref_d, 
                        const EllipticProblem<TypeEquation>& vars, const ElementReference<Dimension, 1>&);
    
    template<class Matrix1, class TypeEquation, class T0>
    static void GetPenalDG(Matrix1& Nabc, const R2& normale, int iquad, int k,
			   int nf, const GlobalGenericMatrix<T0>& nat_mat, int ref, int ref2,
                           const EllipticProblem<TypeEquation>& vars,
			   const ElementReference<Dimension, 1>& Fb);
    
    template<class Vector1, class Vector2, class GenericPb, class T0>
    static void MltPenalDG(const R2& normale, const Vector1& Vn, Vector2& Un,
			   int i, int k, int nf, const GlobalGenericMatrix<T0>& nat_mat,
                           int ref, int ref2, const GenericPb& vars,
			   const ElementReference<Dimension, 1>& Fb);

  };
  
  
  //! time-harmonic Galbrun equation in axisymmetric domain
  class HarmonicGalbrunEquationAxi : public GalbrunAxiEquation_Base<Complex_wp>
  {};


  //! steady Galbrun equation in axisymmetric domain
  class StationaryGalbrunEquationAxi : public GalbrunAxiEquation_Base<Real_wp>
  {};
  

  template<class TypeEquation>  
  class VarGalbrunAxi_Eq : public VarHarmonic<TypeEquation>,
			   public VarGalbrun_Axi<typename TypeEquation::Complexe>
  {
    typedef typename TypeEquation::Complexe Complexe;
    
  public:
    VarGalbrunAxi_Eq();
    
    const R3& GetPhaseOrigin() const;
    bool UseNumericalIntegration(int i) const;
    
    void AddDiracSource(const Real_wp& alpha, Vector<Vector<Real_wp> > & b_source,
			Vector<VirtualSourceFEM<Real_wp, Dimension2>* >& f) const;

    void AddDiracSource(const Complex_wp& alpha, Vector<Vector<Complex_wp> > & b_source,
			Vector<VirtualSourceFEM<Complex_wp, Dimension2>* >& f) const;

    void ModifyOutputUnknown(VectReal_wp& val_u, VectReal_wp& grad_u,
			     int i, const GridInterpolation<Dimension2>& var_interp,
			     int iquad, bool compute_grad) const;
    
    void ModifyOutputUnknown(VectComplex_wp& val_u, VectComplex_wp& grad_u,
			     int i, const GridInterpolation<Dimension2>& var_interp,
			     int iquad, bool compute_grad) const;
    
    void ModifyOutputUnknown(Vector<VectReal_wp>&, Vector<VectReal_wp>&,
                             int, bool, bool) const;
    
    void ModifyOutputUnknown(Vector<VectComplex_wp>&, Vector<VectComplex_wp>&,
                             int, bool, bool) const;

    // additionnal input datas
    void SetInputData(const string& description_field, const VectString& parameters);

    // allocation of arrays containing physical indices
    void InitIndices(int n);

    void GetVaryingIndices(Vector<PhysicalVaryingMedia<Dimension2, Complex_wp>* >& rho_complex,
			   Vector<PhysicalVaryingMedia<Dimension2, Real_wp>* >& rho_real, IVect& num_ref,
			   IVect& num_index, IVect& num_component, Vector<bool>& compute_grad,
			   Vector<bool>& compute_hess);

    void ComputePhysicalCoefficients();

    void AllocateMassMatrices();
    
    void ComputeLocalMassMatrix(int i, int N, bool linear_sparse,
				SetPoints<Dimension2>& PointsElem,
                                SetMatrices<Dimension2>& MatricesElem,
                                IVect& OrderFace, const ElementGeomReference<Dimension2>& Fb);
    
    // set physical indices of media with reference i
    void SetIndices(int i, const VectString& parameters);
    void SetPhysicalIndex(const string& name_media, int i, const VectString& parameters);
    string GetPhysicalIndexName(int m) const;
    int GetNbPhysicalIndices() const;
    bool IsVaryingMedia(int i) const;
    bool IsVaryingMedia(int m, int i) const;
    Real_wp GetVelocityOfMedia(int ref) const;
    Real_wp GetVelocityOfInfinity() const;

    void PerformOtherInitializations();
    void CheckInputMesh();

    void ConstructAll(const string& name_file, const string& name_element,
		      Vector<string>& lines_data_file);
    
  };
  
  
  //! volumetric source
  template<>
  class VolumetricSource<HarmonicGalbrunEquationAxi>
    : public VolumetricSource_AxiGalbrun
  {
  public :
    template<class TypeEquation>
    VolumetricSource(const EllipticProblem<TypeEquation>&,
		     const Vector<VectString>& param);
    
  };

  
  //! class used to solve harmonic Galbrun equation with DG method
  template<>
  class EllipticProblem<HarmonicGalbrunEquationAxi>
    : public VarGalbrunAxi_Eq<HarmonicGalbrunEquationAxi>
  {
  public:    
    void AddElementaryFluxesDG(VirtualMatrix<Real_wp>& mat_sp,
			       const GlobalGenericMatrix<Real_wp>& nat_mat,
			       int offset_row = 0, int offset_col = 0);

    void AddElementaryFluxesDG(VirtualMatrix<Complex_wp>& mat_sp,
			       const GlobalGenericMatrix<Complex_wp>& nat_mat,
			       int offset_row = 0, int offset_col = 0);
    
    void ComputeElementaryMatrix(int, IVect&, VirtualMatrix<Real_wp>&,
				 CondensationBlockSolver_Base<Real_wp>&,
				 const GlobalGenericMatrix<Real_wp>&);

    void ComputeElementaryMatrix(int, IVect&, VirtualMatrix<Complex_wp>&,
				 CondensationBlockSolver_Base<Complex_wp>&,
				 const GlobalGenericMatrix<Complex_wp>&);

  };
  
  
  //! class used to solve real Galbrun equation with DG method
  template<>
  class EllipticProblem<StationaryGalbrunEquationAxi>
    : public VarGalbrunAxi_Eq<StationaryGalbrunEquationAxi>
  {
  public:
    
  };
  
}

#define MONTJOIE_FILE_AXISYM_GALBRUN_HXX
#endif
