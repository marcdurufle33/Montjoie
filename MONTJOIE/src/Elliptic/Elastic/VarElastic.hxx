#ifndef MONTJOIE_FILE_VAR_ELASTIC_HXX

namespace Montjoie
{

  //! class to specify elastic equation for H1 finite element
  template<class T, class Dim>
  class ElasticEquation_Base : public GenericEquation<T>
  {
  public :
    typedef Dim Dimension;
    
    enum {nb_unknowns = Dimension::dim_N, nb_unknowns_scal = Dimension::dim_N,
	  nb_components_en = Dimension::dim_N, nb_components_hn = Dimension::dim_N,
          nb_unknowns_vec = Dimension::dim_N*(Dimension::dim_N+1)/2, nb_unknowns_hdg=Dimension::dim_N};
    
    static inline bool SymmetricGlobalMatrix() { return true; }
    static inline bool SymmetricElementaryMatrix() { return true; }

    template<class TypeEquation>
    static void ComputeMassMatrix(const EllipticProblem<TypeEquation>& var,
                                  int i, const ElementReference_Dim<Dimension>& Fb);
    
    template<class TypeEquation, class T0, class Vector1>
    static void ApplyTensorStiffness(const EllipticProblem<TypeEquation>& var,
                                     int i, int j,
				     const GlobalGenericMatrix<T0>& nat_mat, int ref,
                                     Vector1& dU, Vector1& dV);

    template<class TypeEquation, class T0, class MatStiff>
    static void GetGradGradTensor(const EllipticProblem<TypeEquation>& vars,
				  int num_elem, int jloc, const GlobalGenericMatrix<T0>& nat_mat,
                                  int ref, MatStiff& Cgrad_grad);
    
    template<class Matrix1, class T0, class GenericPb>
    static void GetPenalDG(Matrix1& Nabc, typename Dimension::R_N& normale,
                           int iquad, int k, int nf, const GlobalGenericMatrix<T0>& nat_mat,
			   int ref, int ref2,
                           const GenericPb& vars, const ElementReference<Dim, 1>& Fb);
    
  };
  

  //! class to specify stationary elastic equation for H1 finite element
  template<class Dimension>
  class ElasticEquation : public ElasticEquation_Base<Real_wp, Dimension>
  {
  public :
    
    template<class TypeEquation, class T0, class MatMass>
    static void GetTensorMass(const EllipticProblem<TypeEquation>& vars,
			      int i, int j, const GlobalGenericMatrix<T0>& nat_mat, int ref, MatMass& mass);
    
    template<class TypeEquation, class T0, class Vector1>
    static void ApplyTensorMass(const EllipticProblem<TypeEquation>& var, int i, int j,
				const GlobalGenericMatrix<T0>& nat_mat, int ref, Vector1& U, Vector1& V);
    
  };

  
  //! class to specify time-harmonic elastic equation for H1 finite element
  template<class Dimension>
  class HarmonicElasticEquation : public ElasticEquation_Base<Complex_wp, Dimension>
  {
  public :
    template<class TypeEquation, class T0, class MatMass>
    static void GetTensorMass(const EllipticProblem<TypeEquation>& vars,
			      int i, int j, const GlobalGenericMatrix<T0>& nat_mat, int ref, MatMass& mass);
    
    template<class TypeEquation, class T0, class Vector1>
    static void ApplyTensorMass(const EllipticProblem<TypeEquation>& var, int i, int j,
				const GlobalGenericMatrix<T0>& nat_mat, int ref, Vector1& U, Vector1& V);
    
  };
  
  
  //! base class for solving elastodynamics equations
  class VarElastic_Base
  {
  protected:
    bool display_stress;
    //! coefficient to use in the penalty term
    Vector<Real_wp> coefficient_mu_penalty;
    //! adimensionalization coefficients
    Real_wp omega_bar, u_bar;
    
  private:
    VarProblem_Base& var_problem_base;
    
  public:
    template<class TypeEquation>
    VarElastic_Base(EllipticProblem<TypeEquation>& var);

    void SetInputData(const string& description_field, const VectString& parameters);
    
    template<class Complexe>
    void UpdateShiftAdimensionalization(Complexe&, Complexe&);
    
    template<class Complexe>
    void UpdateEigenvaluesAdimensionalization(Vector<Complexe>&, Vector<Complexe>&, 
					      Matrix<Complexe, General, ColMajor>&);

    template<class T0, class Complexe>
    void ModifyOutputUnknown(const Vector<ElasticPhysicalIndice<Dimension2, 2, T0> >& ref_tensorC,
			     const VarProblem<Dimension2>& var_problem,
			     Vector<Complexe>& val_u,
			     Vector<Complexe>& grad_u, int i,
			     const GridInterpolation<Dimension2>& var_interp,
			     int iquad, bool compute_grad) const;
    
    template<class Complexe>
    void ModifyOutputUnknown(const Vector<ElasticPhysicalIndice<Dimension2, 2, Complexe> >& ref_tensorC,
			     const VarProblem<Dimension2>& var_problem,
			     Vector<Complexe>& val_u,
			     Vector<Complexe>& grad_u, int i,
			     const GridInterpolation<Dimension2>& var_interp,
			     int iquad, bool compute_grad) const;

#ifdef MONTJOIE_WITH_THREE_DIM
    template<class T0, class Complexe>
    void ModifyOutputUnknown(const Vector<ElasticPhysicalIndice<Dimension3, 3, T0> >& ref_tensorC,
			     const VarProblem<Dimension3>& var_problem,
			     Vector<Complexe>& val_u,
			     Vector<Complexe>& grad_u, int i,
			     const GridInterpolation<Dimension3>& var_interp,
			     int iquad, bool compute_grad) const;

    template<class Complexe>
    void ModifyOutputUnknown(const Vector<ElasticPhysicalIndice<Dimension3, 3, Complexe> >& ref_tensorC,
			     const VarProblem<Dimension3>& var_problem,
			     Vector<Complexe>& val_u,
			     Vector<Complexe>& grad_u, int i,
			     const GridInterpolation<Dimension3>& var_interp,
			     int iquad, bool compute_grad) const;
#endif

    template<class Complexe>
    void ComputeEnHnQuadrature(Vector<Vector<Complexe> >& u_quadrature,
			       Vector<Vector<Complexe> >& grad_quadrature,
			       int num_elem, const Vector<R2>& pts, const Vector<R2>& normale,
			       bool compute_H, Vector<Complexe>& En_quad, Vector<Complexe>& Hn_quad) const;

#ifdef MONTJOIE_WITH_THREE_DIM
    template<class Complexe>
    void ComputeEnHnQuadrature(Vector<Vector<Complexe> >& u_quadrature,
			       Vector<Vector<Complexe> >& grad_quadrature,
			       int num_elem, const Vector<R3>& pts, const Vector<R3>& normale,
			       bool compute_H, Vector<Complexe>& En_quad, Vector<Complexe>& Hn_quad) const;
#endif

    template<class Complexe>
    void ComputeEnHnNodal(Vector<Vector<Complexe> >& u_nodal, Vector<Vector<Complexe> >& grad_nodal,
			  int num_elem, const Vector<R2>& pts, const Vector<R2>& normale,
			  Vector<Vector<Complexe> >& En_nodal, Vector<Vector<Complexe> >& Hn_nodal) const;


#ifdef MONTJOIE_WITH_THREE_DIM
    template<class Complexe>
    void ComputeEnHnNodal(Vector<Vector<Complexe> >& u_nodal, Vector<Vector<Complexe> >& grad_nodal,
			  int num_elem, const Vector<R3>& pts, const Vector<R3>& normale,
			  Vector<Vector<Complexe> >& En_nodal, Vector<Vector<Complexe> >& Hn_nodal) const;
#endif
    
    Real_wp GetCoefficientPenaltyStiffness(int ref) const;
    
    template<class T0, class Complexe>
    void ComputeElementaryMatrix(const Vector<ElasticPhysicalIndice<Dimension2, 2, T0> >& ref_tensorC,
				 const Vector<ScalarPhysicalIndice<Dimension2, T0> >& ref_rho,
				 const Vector<ScalarPhysicalIndice<Dimension2, T0> >& ref_sigma,
				 int iquad, IVect& num_dof, VirtualMatrix<Complexe>& mat_interac,
				 const GlobalGenericMatrix<Complexe>& nat_mat,
				 const VarProblem<Dimension2>& var_problem,
				 const VarBoundaryCondition_Dim<T0, Dimension2>& var_boundary,
				 const ElementReference<Dimension2, 1>& Fb);

#ifdef MONTJOIE_WITH_THREE_DIM
    template<class T0, class Complexe>    
    void ComputeElementaryMatrix(const Vector<ElasticPhysicalIndice<Dimension3, 3, T0> >& ref_tensorC,
				 const Vector<ScalarPhysicalIndice<Dimension3, T0> >& ref_rho,
				 const Vector<ScalarPhysicalIndice<Dimension3, T0> >& ref_sigma,
				 int iquad, IVect& num_dof, VirtualMatrix<Complexe>& mat_interac,
				 const GlobalGenericMatrix<Complexe>& nat_mat,
				 const VarProblem<Dimension3>& var_problem,
				 const VarBoundaryCondition_Dim<T0, Dimension3>& var_boundary,
				 const ElementReference<Dimension3, 1>& Fb);
#endif
    
  };
  

  //! class to solve elastics equations
  template<class Complexe, class Dimension>
  class VarElastic_Dim : public VarElastic_Base
  {
  public :
    typedef typename Dimension::R_N R_N; //!< point in R2 or R3
    typedef typename Dimension::VectR_N VectR_N; //!< vector of points in R2 or R3
    typedef Vector<Complexe> VectComplexe; //!< vector of floats or complex floats
    
    //! index rho
    Vector<ScalarPhysicalIndice<Dimension, Complexe> > ref_rho;
    //! damping coefficient
    Vector<ScalarPhysicalIndice<Dimension, Complexe> > ref_sigma;
    //! elastic tensor C
    Vector<ElasticPhysicalIndice<Dimension, Dimension::dim_N, Complexe> > ref_tensorC;
    //! elastic tensor A (inverse of C)
    Vector<ElasticPhysicalIndice<Dimension, Dimension::dim_N, Complexe> > ref_tensorA;

  private:
    VarProblem<Dimension>& var_problem;
    
  public:
    template<class TypeEquation>
    VarElastic_Dim(EllipticProblem<TypeEquation>&);
        
    void InitIndices(int n);
    int GetNbPhysicalIndices() const;
    
    void CopyInputData(const VarElastic_Dim<Complexe, Dimension>& var);

    void SetIndices(int i, const VectString& parameters);
    void SetPhysicalIndex(const string& name_media, int i, const VectString& parameters);
    string GetPhysicalIndexName(int m) const;
    bool IsVaryingMedia(int i) const;
    Real_wp GetVelocityOfMedia(int ref) const;
    Real_wp GetVelocityOfInfinity() const;

    void AllocateMassMatrices();
    
    void GetVaryingIndices(Vector<PhysicalVaryingMedia<Dimension, Complexe>* >& rho_complex,
			   IVect& num_ref, IVect& num_index, IVect& num_component,
			   Vector<bool>& compute_grad, Vector<bool>& compute_hess);
    
    void FinalizeComputationVaryingIndices();
    
    void ComputeNumberOfDofs(int& nodl, IVect& offset_dof_unknown,
			     IVect& OffsetDofV);
    
    void PutOtherGlobalDofs(); 
    
  };


  //! class for solving elastodynamics equation depending on the equation
  template<class TypeEquation>
  class VarElastic_Eq : public VarElastic_Dim<typename TypeEquation::Complexe,
					      typename TypeEquation::Dimension>,
			public VarHarmonic<TypeEquation>
  {
    typedef typename TypeEquation::Dimension Dimension;
    typedef typename TypeEquation::Complexe Complexe;
    typedef typename Dimension::R_N R_N;

  public:    
    VarElastic_Eq();

    void SetInputData(const string& description_field, const VectString& parameters);

    void SetFirstOrderFormulation(bool mix_);
    bool IsSymmetricProblem(bool eigen = false) const;

    bool IsSparseElementaryMatrix(const GlobalGenericMatrix<Real_wp>&) const;
    bool IsSparseElementaryMatrix(const GlobalGenericMatrix<Complex_wp>&) const;

    bool IsSymmetricElementaryMatrix(const GlobalGenericMatrix<Real_wp>&) const;
    bool IsSymmetricElementaryMatrix(const GlobalGenericMatrix<Complex_wp>&) const;

    bool IsDiagonalElementaryMatrix(const GlobalGenericMatrix<Real_wp>&) const;
    bool IsDiagonalElementaryMatrix(const GlobalGenericMatrix<Complex_wp>&) const;
    
    void ComputeNumberOfDofs();
    void PutOtherGlobalDofs(); 
    
    void TreatDirichletCondition();

    void UpdateShiftAdimensionalization(Real_wp&, Real_wp&);    
    void UpdateEigenvaluesAdimensionalization(Vector<Real_wp>&, Vector<Real_wp>&, 
					      Matrix<Real_wp, General, ColMajor>&);

    void UpdateShiftAdimensionalization(Complex_wp&, Complex_wp&);    
    void UpdateEigenvaluesAdimensionalization(Vector<Complex_wp>&, Vector<Complex_wp>&, 
					      Matrix<Complex_wp, General, ColMajor>&);
    
    void ModifyOutputUnknown(Vector<Real_wp>& val_u,
			     Vector<Real_wp>& grad_u, int i,
			     const GridInterpolation<Dimension>& var_interp,
			     int iquad, bool compute_grad) const;

    void ModifyOutputUnknown(Vector<Complex_wp>& val_u,
			     Vector<Complex_wp>& grad_u, int i,
			     const GridInterpolation<Dimension>& var_interp,
			     int iquad, bool compute_grad) const;
    
    void ModifyOutputUnknown(Vector<VectReal_wp>&, Vector<VectReal_wp>&,
                             int, bool, bool) const;
    
    void ModifyOutputUnknown(Vector<VectComplex_wp>&, Vector<VectComplex_wp>&,
                             int, bool, bool) const;

    void ComputeEnHnQuadrature(Vector<VectReal_wp>& u_quadrature,
			       Vector<VectReal_wp>& grad_quadrature,
			       int num_elem, const Vector<R_N>& pts, const Vector<R_N>& normale,
			       bool compute_H, Vector<Real_wp>& En_quad, Vector<Real_wp>& Hn_quad) const;

    void ComputeEnHnQuadrature(Vector<VectComplex_wp>& u_quadrature,
			       Vector<VectComplex_wp>& grad_quadrature,
			       int num_elem, const Vector<R_N>& pts, const Vector<R_N>& normale,
			       bool compute_H, Vector<Complex_wp>& En_quad, Vector<Complex_wp>& Hn_quad) const;

    void ComputeEnHnNodal(Vector<Vector<Real_wp> >& u_nodal, Vector<Vector<Real_wp> >& grad_nodal,
			  int num_elem, const Vector<R_N>& pts, const Vector<R_N>& normale,
			  Vector<Vector<Real_wp> >& En_nodal, Vector<Vector<Real_wp> >& Hn_nodal) const;

    void ComputeEnHnNodal(Vector<Vector<Complex_wp> >& u_nodal, Vector<Vector<Complex_wp> >& grad_nodal,
			  int num_elem, const Vector<R_N>& pts, const Vector<R_N>& normale,
			  Vector<Vector<Complex_wp> >& En_nodal, Vector<Vector<Complex_wp> >& Hn_nodal) const;
    
    Real_wp GetCoefficientPenaltyStiffness(int ref) const;
    
    
    void InitIndices(int n);
    int GetNbPhysicalIndices() const;
    
    void CopyInputData(const VarProblem_Base& var);

    void SetIndices(int i, const VectString& parameters);
    void SetPhysicalIndex(const string& name_media, int i, const VectString& parameters);
    string GetPhysicalIndexName(int m) const;
    bool IsVaryingMedia(int i) const;
    bool IsVaryingMedia(int m, int i) const;
    Real_wp GetVelocityOfMedia(int ref) const;
    Real_wp GetVelocityOfInfinity() const;

    void AllocateMassMatrices();    
    void FinalizeComputationVaryingIndices();
    
  };  
    
  
  //! class for solving elastodynamics equation in time-harmonic domain
  template<class TypeEquation>
  class VarElastic_Cplx : public VarElastic_Eq<TypeEquation>
  {
    typedef typename TypeEquation::Complexe Complexe;
    typedef typename TypeEquation::Dimension Dimension;
    
  public:
    void GetVaryingIndices(Vector<PhysicalVaryingMedia<Dimension, Complex_wp>* >& rho_complex,
			   Vector<PhysicalVaryingMedia<Dimension, Real_wp>* >& rho_real, IVect& num_ref,
                           IVect& num_index, IVect& num_component, Vector<bool>& compute_grad,
			   Vector<bool>& compute_hess);
    
  };


  //! class for solving elastodynamics equation in time-domain
  template<class TypeEquation>
  class VarElastic_Real : public VarElastic_Eq<TypeEquation>
  {
    typedef typename TypeEquation::Complexe Complexe;
    typedef typename TypeEquation::Dimension Dimension;
    
  public:
    void GetVaryingIndices(Vector<PhysicalVaryingMedia<Dimension, Complex_wp>* >& rho_complex,
			   Vector<PhysicalVaryingMedia<Dimension, Real_wp>* >& rho_real, IVect& num_ref,
                           IVect& num_index, IVect& num_component, Vector<bool>& compute_grad,
			   Vector<bool>& compute_hess);
    
  };


  //! class for first-order absorbing boundary condition for elastodynamic equation
  template<class Dimension>
  class ImpedanceABC<Complex_wp, HarmonicElasticEquation<Dimension> >
    : public ImpedanceFunction_Base<Complex_wp, Dimension>
  {
  protected :
    enum{d = Dimension::dim_N};
    
    Vector<TinyMatrix<Complex_wp, Symmetric, d, d> > feval;
    
  private:
    const VarProblem<Dimension>& var_problem;
    const VarElastic_Dim<Complex_wp, Dimension>& var_elas;
    
  public:
    template<class TypeEquation>
    ImpedanceABC(const EllipticProblem<TypeEquation>& var);

    void EvaluateImpedancePhi(int i, int num_elem, int num_edge, int num_loc, int k,
                              const GlobalGenericMatrix<Complex_wp>& nat_mat, int ref,
                              const SetPoints<Dimension>& Pts, const SetMatrices<Dimension>& Mat);
    
    void ApplyImpedancePhi_H1(int n, int j, int, const TinyVector<Real_wp, 1>& phi,
                              const typename Dimension::R_N& dphi, Vector<Complex_wp>& feval_phi);
    
  };


  //! class for first-order absorbing boundary condition for elastodynamic equation
  template<class T, class Dimension>
  class ImpedanceABC<T, ElasticEquation<Dimension> >
    : public ImpedanceFunction_Base<T, Dimension>
  {
  protected :
    enum{d = Dimension::dim_N};
    
    Vector<TinyMatrix<T, Symmetric, d, d> > feval;
    
  private:
    const VarProblem<Dimension>& var_problem;
    const VarElastic_Dim<Real_wp, Dimension>& var_elas;
    
  public:
    template<class TypeEquation>
    ImpedanceABC(const EllipticProblem<TypeEquation>& var);

    void EvaluateImpedancePhi(int i, int num_elem, int num_edge, int num_loc, int k,
                              const GlobalGenericMatrix<T>& nat_mat, int ref,
                              const SetPoints<Dimension>& Pts, const SetMatrices<Dimension>& Mat);
    
    void ApplyImpedancePhi_H1(int n, int j, int, const TinyVector<Real_wp, 1>& phi,
                              const typename Dimension::R_N& dphi, Vector<T>& feval_phi);
    
  };

  
  //! class used to solve real elastic equation with H1 elements
  template<class Dimension>
  class EllipticProblem<ElasticEquation<Dimension> >
    : public VarElastic_Real<ElasticEquation<Dimension> >
  {
  public:
    void ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Real_wp>& mat_elem,
				 CondensationBlockSolver_Base<Real_wp>&,
				 const GlobalGenericMatrix<Real_wp>& nat_mat);

    void ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Complex_wp>& mat_elem,
				 CondensationBlockSolver_Base<Complex_wp>&,
				 const GlobalGenericMatrix<Complex_wp>& nat_mat);

    void AddElementaryFluxesDG(VirtualMatrix<Real_wp>& mat_sp,
			       const GlobalGenericMatrix<Real_wp>& nat_mat,
			       int offset_row = 0, int offset_col = 0);

    void AddElementaryFluxesDG(VirtualMatrix<Complex_wp>& mat_sp,
			       const GlobalGenericMatrix<Complex_wp>& nat_mat,
			       int offset_row = 0, int offset_col = 0);
    
  };


  template<class Dimension>
  class FemMatrixFreeClass<Real_wp, ElasticEquation<Dimension> >
    : public FemMatrixFreeClass_Eq<Real_wp, ElasticEquation<Dimension> >
  {
  public:
    
    FemMatrixFreeClass(const EllipticProblem<ElasticEquation<Dimension> >& var_);
    
    virtual void MltAddFree(const GlobalGenericMatrix<Real_wp>& nat_mat,
			    const SeldonTranspose&, int lvl, 
			    const Vector<Real_wp>& X, Vector<Real_wp>& Y) const;

    virtual void MltAddFree(const GlobalGenericMatrix<Real_wp>& nat_mat,
			    const SeldonTranspose&, int lvl, 
			    const Vector<Complex_wp>& X, Vector<Complex_wp>& Y) const;
    
  };


  //! class used to solve time-harmonic elastic equation with H1 elements
  template<class Dimension>
  class EllipticProblem<HarmonicElasticEquation<Dimension> >
    : public VarElastic_Cplx<HarmonicElasticEquation<Dimension> >
  {
  public:
    void ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Real_wp>& mat_elem,
				 CondensationBlockSolver_Base<Real_wp>&,
				 const GlobalGenericMatrix<Real_wp>& nat_mat);

    void ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Complex_wp>& mat_elem,
				 CondensationBlockSolver_Base<Complex_wp>&,
				 const GlobalGenericMatrix<Complex_wp>& nat_mat);
    
    void AddElementaryFluxesDG(VirtualMatrix<Real_wp>& mat_sp,
			       const GlobalGenericMatrix<Real_wp>& nat_mat,
			       int offset_row = 0, int offset_col = 0);

    void AddElementaryFluxesDG(VirtualMatrix<Complex_wp>& mat_sp,
			       const GlobalGenericMatrix<Complex_wp>& nat_mat,
			       int offset_row = 0, int offset_col = 0);
    
  };


  template<class Dimension>
  class FemMatrixFreeClass<Complex_wp, HarmonicElasticEquation<Dimension> >
    : public FemMatrixFreeClass_Eq<Complex_wp, HarmonicElasticEquation<Dimension> >
  {
  public:
    
    FemMatrixFreeClass(const EllipticProblem<HarmonicElasticEquation<Dimension> >& var_);

    virtual void MltAddFree(const GlobalGenericMatrix<Complex_wp>& nat_mat,
			    const SeldonTranspose&, int lvl, 
			    const Vector<Real_wp>& X, Vector<Real_wp>& Y) const;
    
    virtual void MltAddFree(const GlobalGenericMatrix<Complex_wp>& nat_mat,
			    const SeldonTranspose&, int lvl, 
			    const Vector<Complex_wp>& X, Vector<Complex_wp>& Y) const;
    
  };
  

  //! Elastic equation with Local Discontinuous Galerkin 
  /*!
    In 2-D, these equations are given as :
    du_x/dt - d/dx(sigma_xx) - d/dy(sigma_xy) = 0
    du_y/dt - d/dx(sigma_xy) - d/dy(sigma_yy) = 0
    dsigma/dt - C nabla u = 0
  */
  template<class T, class Dim>
  class ElasticEquationDG_Base : public GenericEquation<T>
  {
  public :
    typedef Dim Dimension;
    
    static const bool FirstOrderFormulation = true;
    
    enum {nb_unknowns = Dimension::dim_N*(Dimension::dim_N+3)/2,
	  nb_components_en = Dimension::dim_N, nb_components_hn = Dimension::dim_N,
          nb_unknowns_scal = nb_unknowns, nb_unknowns_vec = 0, nb_unknowns_hdg=Dimension::dim_N};

    static inline bool SymmetricGlobalMatrix() { return false; }
    static inline bool SymmetricElementaryMatrix() { return false; }
        
    // for compatbility purpose
    template<class TypeEquation>
    static void ComputeMassMatrix(EllipticProblem<TypeEquation>& var,
                                  int num_elem, const ElementReference_Dim<Dimension>& Fb);
        
    // For a detailed description of the following methods
    // look at the file GenericEquation.hxx (class GenericEquation_Base)
    template<class GenericPb, class T0, class Matrix1>
    static void GetNeededDerivative(const GenericPb& vars,
                                    const GlobalGenericMatrix<T0>& nat_mat, Matrix1& At, Matrix1& A);
        
    template<class TypeEquation, class T0, class MatStiff>
    static void GetGradPhiTensor(const EllipticProblem<TypeEquation>& vars,
                                 int num_elem, int jloc, const GlobalGenericMatrix<T0>& nat_mat, int ref,
                                 MatStiff& Dgrad_phi, MatStiff& Ephi_grad);
    
    template<class TypeEquation, class T0>
    static void ApplyGradientUnknown(const EllipticProblem<TypeEquation>& var,
				     int i, int j, const GlobalGenericMatrix<T0>& nat_mat,
                                     int ref, const TinyVector<TinyVector<T0, 2>, 5>& dU,
                                     TinyVector<T0, 5>& V);

    template<class TypeEquation, class T0>
    static void ApplyGradientUnknown(const EllipticProblem<TypeEquation>& var,
				     int i, int j, const GlobalGenericMatrix<T0>& nat_mat,
                                     int ref, const TinyVector<TinyVector<T0, 3>, 9>& dU,
                                     TinyVector<T0, 9>& V);

    template<class TypeEquation, class T0>
    static void ApplyGradientFctTest(const EllipticProblem<TypeEquation>& var,
				     int i, int j, const GlobalGenericMatrix<T0>& nat_mat,
                                     int ref, const TinyVector<T0, 5>& U,
                                     TinyVector<TinyVector<T0, 2>, 5>& dV);

    template<class TypeEquation, class T0>
    static void ApplyGradientFctTest(const EllipticProblem<TypeEquation>& var,
				     int i, int j, const GlobalGenericMatrix<T0>& nat_mat,
                                     int ref, const TinyVector<T0, 9>& U,
                                     TinyVector<TinyVector<T0, 3>, 9>& dV);
    
    template<class TypeEquation, class T0, class MatMass>
    static void GetTensorMass(const EllipticProblem<TypeEquation>& vars,
			      int i, int j, const GlobalGenericMatrix<T0>& nat_mat, int ref, MatMass& mass);
    
    template<class TypeEquation, class T0>
    static void ApplyTensorMass(const EllipticProblem<TypeEquation>& var, int i, int j,
				const GlobalGenericMatrix<T0>& nat_mat, int ref,
                                const TinyVector<T0, 5>& U, TinyVector<T0, 5>& V);

    template<class TypeEquation, class T0>
    static void ApplyTensorMass(const EllipticProblem<TypeEquation>& var, int i, int j,
				const GlobalGenericMatrix<T0>& nat_mat, int ref,
                                const TinyVector<T0, 9>& U, TinyVector<T0, 9>& V);

    template<class Matrix1, class TypeEquation, class T0>
    static void GetNabc(Matrix1& Nabc, const typename Dimension::R_N& normale,
                        int ref, int num_elem, int jloc,
                        const GlobalGenericMatrix<T0>& nat_mat, int ref_d,
                        const EllipticProblem<TypeEquation>& vars, const ElementReference<Dimension, 1>&);
    
    template<class Vector1, class TypeEquation, class T0>
    static void MltNabc(typename Dimension::R_N& normale, int ref, const Vector1& Vn, Vector1& Un,
                        int num_elem1, int num_point, const GlobalGenericMatrix<T0>& nat_mat, int ref_d,
                        const EllipticProblem<TypeEquation>& vars, const ElementReference<Dimension, 1>& Fb);

    template<class Matrix1, class GenericPb, class T0>
    static void GetPenalDG(Matrix1& Nabc, typename Dimension::R_N& normale, int iquad, int k,
                           int nf, const GlobalGenericMatrix<T0>& nat_mat, int ref, int ref2,
                           const GenericPb& vars, const ElementReference<Dimension, 1>& Fb);
    
    template<class T0, class GenericPb>
    static void MltPenalDG(const typename Dimension::R_N& normale,
                           const TinyVector<T0, 5>& Vn, TinyVector<T0, 5>& Un,
                           int i, int k, int nf, const GlobalGenericMatrix<T0>& nat_mat,
                           int ref, int ref2, const GenericPb& vars,
			   const ElementReference<Dimension, 1>& Fb);

    template<class T0, class GenericPb>
    static void MltPenalDG(const typename Dimension::R_N& normale,
                           const TinyVector<T0, 9>& Vn, TinyVector<T0, 9>& Un,
                           int i, int k, int nf, const GlobalGenericMatrix<T0>& nat_mat,
                           int ref, int ref2, const GenericPb& vars,
			   const ElementReference<Dimension, 1>& Fb);
    
  };
  
  
  //! static elastic equation with LDG
  template<class Dim>
  class ElasticEquationDG
    : public ElasticEquationDG_Base<Real_wp, Dim>
  {
  };


  //! time-harmonic elastic equation with LDG
  template<class Dim>
  class HarmonicElasticEquationDG
    : public ElasticEquationDG_Base<Complex_wp, Dim>
  {
  };
  

  //! class used to solve real elastic equation with DG
  template<class Dimension>
  class EllipticProblem<ElasticEquationDG<Dimension> >
    : public VarElastic_Real<ElasticEquationDG<Dimension> >
  {
  public:
    void ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Real_wp>& mat_elem,
				 CondensationBlockSolver_Base<Real_wp>&,
				 const GlobalGenericMatrix<Real_wp>& nat_mat);

    void ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Complex_wp>& mat_elem,
				 CondensationBlockSolver_Base<Complex_wp>&,
				 const GlobalGenericMatrix<Complex_wp>& nat_mat);

    void AddElementaryFluxesDG(VirtualMatrix<Real_wp>& mat_sp,
			       const GlobalGenericMatrix<Real_wp>& nat_mat,
			       int offset_row = 0, int offset_col = 0);

    void AddElementaryFluxesDG(VirtualMatrix<Complex_wp>& mat_sp,
			       const GlobalGenericMatrix<Complex_wp>& nat_mat,
			       int offset_row = 0, int offset_col = 0);
  };
  
  
  //! class used to solve time-harmonic elastic equation with DG elements
  template<class Dimension>
  class EllipticProblem<HarmonicElasticEquationDG<Dimension> >
    : public VarElastic_Cplx<HarmonicElasticEquationDG<Dimension> >
  {
  public:
    void ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Real_wp>& mat_elem,
				 CondensationBlockSolver_Base<Real_wp>&,
				 const GlobalGenericMatrix<Real_wp>& nat_mat);

    void ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Complex_wp>& mat_elem,
				 CondensationBlockSolver_Base<Complex_wp>&,
				 const GlobalGenericMatrix<Complex_wp>& nat_mat);
    
    void AddElementaryFluxesDG(VirtualMatrix<Real_wp>& mat_sp,
			       const GlobalGenericMatrix<Real_wp>& nat_mat,
			       int offset_row = 0, int offset_col = 0);

    void AddElementaryFluxesDG(VirtualMatrix<Complex_wp>& mat_sp,
			       const GlobalGenericMatrix<Complex_wp>& nat_mat,
			       int offset_row = 0, int offset_col = 0);
    
  };
  
}

#define MONTJOIE_FILE_VAR_ELASTIC_HXX
#endif
