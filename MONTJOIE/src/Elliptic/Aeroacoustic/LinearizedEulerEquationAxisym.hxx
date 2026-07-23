#ifndef MONTJOIE_FILE_LINEARIZED_EULER_EQUATION_AXISYM_HXX

namespace Montjoie
{
  
  //! class to specify aeroacoustic equation for any flow in cylindrical coordinates
  /*!
    Considered equations are
    (-i omega + sigma + m_r d/dr -i m m_theta/r + m_z d/dz) p 
          + 1/r d/dr(r c_0^2 u_r) - i m/r c_0^2 u_theta + d/dz(c_0^2 u_z) = 0
    
    (-i omega + sigma + m_r d/dr -i m m_theta/r + m_z d/dz) u_r + dp/dr + dm_r/dr u_r + dm_r/dz u_z = 0
    (-i omega + sigma + m_r d/dr -i m m_theta/r + m_z d/dz) u_theta
                - im/r p + dm_\theta/dr u_r + dm_\theta/dz u_z = 0
    (-i omega + sigma + m_r d/dr -i m m_theta/r + m_z d/dz) u_z + dp/dz + dm_z/dr u_r + dm_z/dz u_z = 0
    
    If Bogey-Bailly-Juve model is selected, the terms with derivatives of m_r and m_z are dropped
    If conservative model is selected, coefficients are modified (see Aeroacoustic.hxx)
    For Galbrun model, we have the following equations:
    rho_0 (-i omega + sigma + m_r d/dr -i m m_theta/r + m_z d/dz) p 
          + \rho_0^2 c_0^2 ( 1/r d/dr(r u_r) - i m/r u_theta + d/dz(u_z) ) = 0
     
    rho_0 (-i omega + sigma + m_r d/dr -i m m_theta/r + m_z d/dz) u_r + dp/dr
             + [ (u_r/r -i m u_theta/r + du_z/dz) dp0/dr - du_z/dr dp0/dz ] / (-i omega + sigma) = 0

    rho_0 (-i omega + sigma + m_r d/dr -i m m_theta/r + m_z d/dz) u_theta - im p/r
             + [ (i m u_r/r dp0/dr + i m u_z/r dp0/dz) ] / (-i omega + sigma) = 0

    rho_0 (-i omega + sigma + m_r d/dr -i m m_theta/r + m_z d/dz) u_z + dp/dz
             + [ (u_r/r -i m u_theta/r + du_r/dr) dp0/dz - du_r/dz dp0/dr ] / (-i omega + sigma) = 0
   */
  template<class T>
  class AeroacousticAxiEquation_Base : public GenericEquation<T>
  {
  public :
    typedef Dimension2 Dimension;
    
    static const bool FirstOrderFormulation = true;
    
    enum {nb_unknowns = 4, nb_components_u = 1,
          nb_components_grad = 2, nb_unknowns_hdg=1,
	  nb_components_en = 1, nb_components_hn = 1,
	  nb_unknowns_scal = 4, nb_unknowns_vec = 0};

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
    static void GetAbsoluteD(TinyMatrix<T0, General, 4, 4>& Dtest, const R2& normale, const Real_wp& alpha,
                             const Real_wp& c0, const T0& rtilde, bool conservative, bool galbrun, const T0& gamma);
    
    template<class Matrix1, class GenericPb, class T0>
    static void GetAbsoluteD(Matrix1& Nabc, const R2& normale, int iquad, int k,
                             const GlobalGenericMatrix<T0>& nat_mat, const GenericPb& vars);
    
    template<class Matrix1, class GenericPb, class T0>
    static void GetAbsoluteMatrixD(Matrix1& D, const R2& normale, int iquad, int k,
                                   const GlobalGenericMatrix<T0>& nat_mat, const GenericPb& vars);      
    
    template<class Matrix1, class GenericPb, class T0>
    static void GetNabc(Matrix1& Nabc, const R2& normale,
			int ref, int iquad, int k, const GlobalGenericMatrix<T0>& nat_mat, int ref_d, 
			const GenericPb& vars, const ElementReference<Dimension, 1>& Fb);
    
    template<class Vector1, class TypeEquation, class T0>
    static void MltNabc(const R2& normale, int ref,
                        const Vector1& Vn, Vector1& Un, int num_elem1,
                        int num_point, const GlobalGenericMatrix<T0>& nat_mat, int ref_d, 
                        const EllipticProblem<TypeEquation>& vars, const ElementReference<Dimension, 1>& );
    
    template<class Matrix1, class TypeEquation, class T0>
    static void GetPenalDG(Matrix1& Nabc, const R2& normale, int iquad, int k,
			   int nf, const GlobalGenericMatrix<T0>& nat_mat, int ref, int ref2,
                           const EllipticProblem<TypeEquation>& vars,
			   const ElementReference<Dimension, 1>& Fb);
    
    template<class Vector1, class Vector2, class GenericPb, class T0>
    static void MltPenalDG(const R2& normale, const Vector1& Vn, Vector2& Un,
			   int i, int k, int nf, const GlobalGenericMatrix<T0>& nat_mat,
                           int ref, int ref2, const GenericPb& vars, const ElementReference<Dimension, 1>& Fb);
    
  };

  
  //! in time-harmonic domain
  class HarmonicAeroacousticAxiEquation : public AeroacousticAxiEquation_Base<Complex_wp>
  {};

  
  //! stationary equation
  class StationaryAeroacousticAxiEquation : public AeroacousticAxiEquation_Base<Real_wp>
  {};

  
  //! Class specifying Linearized Euler Equations (LEE) in cylindrical coordinates
  /*!
  Considered equations are
    (-i omega + sigma + m_r d/dr -i m m_theta/r + m_z d/dz) p 
          + 1/r d/dr(r c_0^2 u_r) - i m/r c_0^2 u_theta + d/dz(c_0^2 u_z) 
          + gamma (div M) p - (gamma-1)/rho_0 (u_r dp0/dr + u_z dp0/dz ) = 0
    
    (-i omega + sigma + m_r d/dr -i m m_theta/r + m_z d/dz) u_r
         + (div M) u_r + dp/dr + dm_r/dr (u_r + rho m_r) + dm_r/dz (u_z + rho m_z) = 0
    
    (-i omega + sigma + m_r d/dr -i m m_theta/r + m_z d/dz) u_theta
          + (div M) u_theta - im/r p + dm_\theta/dr (u_r + rho m_r) + dm_\theta/dz (u_z + rho m_z) = 0
    
    (-i omega + sigma + m_r d/dr -i m m_theta/r + m_z d/dz) u_z 
         + (div M) u_z + dp/dz + dm_z/dr (u_r + rho m_r) + dm_z/dz (u_z + rho m_z) = 0
		
   (-i omega + sigma + m_r d/dr -i m m_theta/r + m_z d/dz) rho
                + rho (div M) + 1/r d/dr(r u_r) - i m/r u_\theta + du_z/dz = 0
  */
  template<class T>
  class LinearizedEulerEquation_Axi : public GenericEquation<T>
  {
  public :
    typedef Dimension2 Dimension;
    
    static const bool FirstOrderFormulation = true;
    
    enum {nb_unknowns = 5, nb_components_u = 1,
          nb_components_grad = 2, nb_unknowns_hdg=0,
	  nb_components_en = 1, nb_components_hn = 1,
	  nb_unknowns_scal = 5, nb_unknowns_vec = 0};

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
    
    template<class Matrix1, class GenericPb, class T0>
    static void GetNabc(Matrix1& Nabc, const R2& normale,
			int ref, int iquad, int k, const GlobalGenericMatrix<T0>& nat_mat, int ref_d, 
			const GenericPb& vars, const ElementReference<Dimension, 1>& Fb);
    
    template<class Vector1, class TypeEquation, class T0>
    static void MltNabc(const R2& normale, int ref,
                        const Vector1& Vn, Vector1& Un, int num_elem1,
                        int num_point, const GlobalGenericMatrix<T0>& nat_mat, int ref_d, 
                        const EllipticProblem<TypeEquation>& vars, const ElementReference<Dimension, 1>& );
    
    template<class Matrix1, class TypeEquation, class T0>
    static void GetPenalDG(Matrix1& Nabc, const R2& normale, int iquad, int k,
			   int nf, const GlobalGenericMatrix<T0>& nat_mat, int ref, int ref2,
                           const EllipticProblem<TypeEquation>& vars, const ElementReference<Dimension, 1>& Fb);
    
    template<class Vector1, class Vector2, class GenericPb, class T0>
    static void MltPenalDG(const R2& normale, const Vector1& Vn, Vector2& Un,
			   int i, int k, int nf, const GlobalGenericMatrix<T0>& nat_mat,
                           int ref, int ref2, const GenericPb& vars,
			   const ElementReference<Dimension, 1>& Fb);
    
  };

  
  //! in time-harmonic domain
  class HarmonicLinearizedEulerEquationAxi : public LinearizedEulerEquation_Axi<Complex_wp>
  {};

  
  //! stationary equation
  class StationaryLinearizedEulerEquationAxi : public LinearizedEulerEquation_Axi<Real_wp>
  {};
  

  template<class Complexe>
  class VarAeroacoustic_Axi : public VarGalbrun_Axi<Complexe>
  {
  private:
    DistributedProblem<Dimension2>& var_problem;
    
  public:
    template<class TypeEquation>
    VarAeroacoustic_Axi(EllipticProblem<TypeEquation>&);

    void SetInputData(const string& description_field, const VectString& parameters);
    
  };
  

  //! class to solve Linearized Euler equations on axisymmetric domains
  template<class TypeEquation>
  class VarAeroacousticAxi_Eq : public VarHarmonic<TypeEquation>,
			      public VarAeroacoustic_Axi<typename TypeEquation::Complexe>
  {
  public :
    typedef typename TypeEquation::Complexe Complexe;
    
    VarAeroacousticAxi_Eq();

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
  class VolumetricSource_AxiAero
    : public VirtualSourceFEM<Complex_wp, Dimension2>
  {
  private:
    //const DistributedProblem<Dimension2>& var_problem;
    const VarAeroacoustic_Axi<Complex_wp>& var_aero;
    const VarBoundaryCondition_Dim<Complex_wp, Dimension2>& var_boundary;

    GaussianSource<Dimension3> fsrc;    
    Vector<Complex_wp> polar;
    
  public :
    Complex_wp coef_vol;
    
    //! constructor with given problem
    template<class TypeEquation>
    VolumetricSource_AxiAero(const EllipticProblem<TypeEquation>& var,
			     const Vector<VectString>& param); 
    
    //! volumetric source
    bool IsNonNullVolumetricSource(const VectR2& s);
       
    //! evaluation of volumetric source
    void EvaluateVolumetricSource(int i, int j, const R2& x, Vector<Complex_wp>& f);
    
  };
  

  //! volumetric source
  template<>
  class VolumetricSource<HarmonicAeroacousticAxiEquation>
    : public VolumetricSource_AxiAero
  {
  public :
    template<class TypeEquation>
    VolumetricSource(const EllipticProblem<TypeEquation>& var,
		     const Vector<VectString>& param);
     
  };


  //! volumetric source
  template<>
  class VolumetricSource<HarmonicLinearizedEulerEquationAxi>
    : public VolumetricSource_AxiAero
  {
  public :
    template<class TypeEquation>
    VolumetricSource(const EllipticProblem<TypeEquation>& var,
		     const Vector<VectString>& param);
     
  };

  
  template<>
  class EllipticProblem<HarmonicAeroacousticAxiEquation>
    : public VarAeroacousticAxi_Eq<HarmonicAeroacousticAxiEquation>
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


  template<>
  class EllipticProblem<StationaryAeroacousticAxiEquation>
    : public VarAeroacousticAxi_Eq<StationaryAeroacousticAxiEquation>
  {
  public:
    
  };


  template<>
  class EllipticProblem<HarmonicLinearizedEulerEquationAxi>
    : public VarAeroacousticAxi_Eq<HarmonicLinearizedEulerEquationAxi>
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


  template<>
  class EllipticProblem<StationaryLinearizedEulerEquationAxi>
    : public VarAeroacousticAxi_Eq<StationaryLinearizedEulerEquationAxi>
  {
  public:
    
  };
  
};

#define MONTJOIE_FILE_LINEARIZED_EULER_EQUATION_AXISYM_HXX
#endif
