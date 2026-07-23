#ifndef MONTJOIE_FILE_VAR_GEOMETRY_PROBLEM_HXX

namespace Montjoie
{
  
  //! class handling mesh and geometry
  template<class Dimension>
  class VarGeometryProblem : public VarProblem_Base
  {
  public :
    typedef typename Dimension::R_N R_N; //!< point in R2 or R3
    typedef typename Dimension::VectR_N VectR_N; //!< vector of points 
    typedef typename Dimension::MatrixN_N MatrixN_N; //!< 2x2 matrix or 3x3 matrix
    typedef typename Dimension::DimensionBoundary DimensionBoundary;

    // Mesh used
    Mesh<Dimension> mesh; //!< mesh related to the problem
    
    //! all quadrature points of the mesh (volume and boundary)
    Vector<VectR_N> Glob_PointsQuadrature;
    bool write_quadrature_points, write_quad_points_pml;
    
    //! jacobian multiplied by weight of integration on each quadrature point of the mesh
    /*!
      Glob_jacobian(i)(j) will be \f$ J_i(\xi_j) \omega_j \f$
      This array is constructed by ComputeMassMatrix() method
    */
    Vector<VectReal_wp> Glob_jacobian;
    //! decomposition of jacobian on polynomial basis
    Vector<VectReal_wp> Glob_decomp_jacobian;

    //! normale at quadrature points of faces/edges
    /*!
      Glob_normale(i,j) will be \f$ n_i(\xi_j) \f$
      This array if constructed by ComputeMassMatrix method and only for dg method
    */
    Vector<VectR_N> Glob_normale;
    //! surface integration element
    Vector<VectReal_wp> Glob_dsj;

    //! element is orthogonal ?
    Vector<int> OrthogonalElement;
    
    //! inverse of jacobian matrices, multiplied by jacobian
    /*!
      Glob_DFjm1(i)(j) will be \f$ J_i(\xi_j) DF_i^{-1}(\xi_j) \f$
      This array is constructed by ComputeMassMatrix() method
    */
    Vector<Vector<MatrixN_N> > Glob_DFjm1;
    //! array to know is the face is "original"
    /*!
      IsNewFace(i,j) is true if it is the first time
      that the matching global face is found
    */
    Vector<Vector<bool> > IsNewFace;

  protected :
    //! wave vector k 
    R_N kwave;
    //! polarization
    Vector<Real_wp> polarization, polarization_grad;
    R_N origine_phase; //!< phase origin in the definition of plane wave
    int reference_infinity; //!< reference of the infinite medium for stratified media

    MeshNumbering<Dimension> mesh_num; //!< numbering of the main unknown
    Vector<MeshNumbering<Dimension>* > other_mesh_num; //!< other numberings (if needed)
    
  private :
    VarBoundaryCondition_Base& var_boundary;
    DistributedProblem_Base& var_comm;
    VarOutputProblem_Dim<Dimension>& var_output;

    void InitDefaultValues();

  public :
    template<class TypeEquation>
    VarGeometryProblem(EllipticProblem<TypeEquation>&);

    ~VarGeometryProblem();
    
    int GetReferenceInfinity() const;
    // Inline methods
    
    const R_N& GetWaveVector() const;
    void SetWaveVector(const R_N&);
    const R_N& GetPhaseOrigin() const;
    const VectReal_wp& GetPolarization() const;
    const VectReal_wp& GetPolarizationGrad() const;
    
    template<class T> void GetPolarization(Vector<T>&) const;
    template<class T> void GetPolarizationGrad(Vector<T>&) const;    
    template<class T, int p> void GetPolarization(TinyVector<T, p>&) const;
    
    template<int p>
    void SetPolarization(const TinyVector<Real_wp, p>& P);    

    virtual void FinalizeComputationVaryingIndices();
    virtual void AllocateMassMatrices();
    
    void DoNotComputeGrid();

    int GetDimension() const;
    
    template<class Vector1, class Vector2>
    void GetLocalUnknownVector(const Vector1& B, int i, Vector2& Uloc) const;

    template<class T0, class T1, int nb_u, class Vector2>
    void AddLocalUnknownVector(const T0& alpha, TinyVector<Vector<T1>, nb_u>& Uloc,
			       int i, Vector2& C) const;

    template<class Vector1>
    void ModifyLocalUnknownVector(Vector1& U, int iquad) const;
    
    template<class Vector1>
    void ModifyLocalComponentVector(Vector1& Eloc, int iquad) const;
    
    template<class Vector1>
    void GetGlobalUnknownVector(Vector1& U, int iquad) const;

    bool InsidePML(int i) const;
    virtual bool ElementInsidePML(int) const;

    virtual void WriteMesh(const string&);
    virtual void SetSameNumberPeriodicDofs();

    // Convenient methods
  protected:
    void UpdateWaveVector();
    void InitPolarization();
    void FillWaveVectorComponents(Real_wp&, Real_wp&, Real_wp&);

  public:
    int GetNbComponentsType(int type) const;
    int GetNbComponentsGradType(int type) const;
    int GetNbComponentsUnknown(int n) const;
    int GetNbComponentsGradient(int n) const;
    int GetNbComponentsAll(int nb_vec = 0) const;
    int GetNbComponentsGradientAll(int nb_vec = 0) const;
    int GetNbComponentsHessianAll(int nb_vec = 0) const;
    
    const MeshNumbering<Dimension>& GetMeshNumbering(int n) const;
    MeshNumbering<Dimension>& GetMeshNumbering(int n);

    void SetInputData(const string& description_field, const VectString& parameters);
    bool FaceHasToBeConsideredForBoundaryIntegral(int num_face) const;
    
    void GetMemoryUsed(map<string, size_t>& var) const;

    void WriteNodalPointsMesh();

    void LocalizePointsBoundaryElement(const VectR_N& Points, int, int,
				       FjInverseProblem<Dimension>& inverseFj,
				       typename DimensionBoundary::VectR_N& pts, IVect& num) const;
    
    bool UseNumericalIntegration(int i) const;

    Real_wp GetWeightedJacobian(int num_elem, int j, bool affine,
                                const ElementGeomReference<Dimension>& Fb) const;

    Real_wp GetSurfaceWeightedJacobian(int num_elem, int num_loc, int num_edge, int j,
                                       const ElementGeomReference<Dimension>& Fb) const;
    
    void GetInverseJacobianMatrix(int num_elem, int j, bool affine,
				  const ElementGeomReference<Dimension>& Fb,
				  MatrixN_N& dfjm1, Real_wp& jacob, Real_wp& jacob_weighted) const;
    
    void FillQuadratureJacobian(int i, int N, const ElementGeomReference<Dimension>& Fb,
				const Vector<R_N>& s, SetPoints<Dimension>& PointsElem,
				SetMatrices<Dimension>& MatricesElem) const;
    
    virtual void ComputeLocalMassMatrix(int i, int nb_points_quadrature, bool linear_sparse_mass,
					SetPoints<Dimension>& PointsElem,
					SetMatrices<Dimension>& MatricesElem,
					IVect& OrderFace, const ElementGeomReference<Dimension>& Fb);
    
    virtual void ComputeMassMatrix(bool compute_rho = true, bool delete_points = true);
    void ClearMassMatrix();
    void ComputeVariableOrder();
    
    Mesh<Dimension>& 
    GetPhysicalCoefficientMesh(int i, const string& file_name, int r, bool same_mesh,
			       Vector<string>& all_names, Vector<Mesh<Dimension> >& all_mesh);

    Vector<FiniteElementInterpolator>& 
    GetPhysicalCoefInterp(int r, const Mesh<Dimension>& meshb,
			  const TinyVector<IVect, 4>& order_mesh,
			  Vector<Vector<FiniteElementInterpolator> >& all_interp,
			  Vector<int> & order_interp);
    
    Vector<FiniteElementInterpolator>&
    GetPhysicalCoefInterpSurf(int r, const Mesh<Dimension>& meshb,
			      const TinyVector<IVect, 4>& order_mesh,
			      Vector<Vector<FiniteElementInterpolator> >& all_interp,
			      Vector<int> & order_interp);

  protected:
    template<class T>
    void ComputePhysicalCoefficientsRef(int ind, int ref, int nb_elt_on_ref,
					int num_index, int num_component,
					bool compute_grad, bool compute_hess,
					Vector<VectR_N>& PointsQuadrature,
					PhysicalVaryingMedia<Dimension, T>& fct_rho,
					Vector<Mesh<Dimension> >& all_mesh,
					Vector<string>& all_names,
					Vector<Vector<FiniteElementInterpolator> >& all_interp,
					Vector<int>& order_interp,
					Vector<Vector<FiniteElementInterpolator> >& all_interp_surf,
					Vector<int>& order_interp_surf, const string& var_file_name);

  public:
    virtual void ComputePhysicalCoefficients();

    virtual void CheckInputMesh();
    virtual void ComputeMeshAndFiniteElement(const string& name_elt, bool split = true);
    
    void CopyInputData(const VarProblem_Base& var);

  protected:
    // Methods to check continuity
    void CheckLocalContinuity(int i, int num_elem, int num_elem2,
			      SetPoints<Dimension>&, SetPoints<Dimension>&,
			      SetMatrices<Dimension>&, SetMatrices<Dimension>&,
			      const ElementReference<Dimension, 3>& Fb,
			      const ElementReference<Dimension, 3>& Fb2);
    
    void CheckLocalContinuity(int i, int num_elem, int num_elem2,
			      SetPoints<Dimension>&, SetPoints<Dimension>&,
			      SetMatrices<Dimension>&, SetMatrices<Dimension>&,
			      const ElementReference<Dimension, 2>& Fb,
			      const ElementReference<Dimension, 2>& Fb2);
        
    void CheckLocalContinuity(int i, int num_elem, int num_elem2,
                              SetPoints<Dimension>&, SetPoints<Dimension>&,
			      SetMatrices<Dimension>&, SetMatrices<Dimension>&,
                              const ElementReference<Dimension, 1>&,
			      const ElementReference<Dimension, 1>&);    

  public:
    // Pure virtual methods
    virtual void GetVaryingIndices(Vector<PhysicalVaryingMedia<Dimension, Complex_wp>* >& rho_complex,
				   Vector<PhysicalVaryingMedia<Dimension, Real_wp>* >& rho_real,
				   IVect& num_ref, IVect& num_index, IVect& num_component,
				   Vector<bool>& compute_grad, Vector<bool>& compute_hess) = 0;    

    virtual void ComputeStoreCoefficientsPML(int i1, int ne, const VectR_N& AllPoints) = 0;
    
    virtual VectR_N PointsQuadInsideND(int i) const = 0;
    virtual const typename DimensionBoundary::VectR_N& PointsQuadratureBoundary(int, int) const = 0;    
    virtual const typename DimensionBoundary::VectR_N& PointsDofBoundary(int, int) const = 0;    
    virtual const ElementGeomReference<Dimension>& GetShapeElement(int i) const = 0;
    
  protected:
    virtual void ComputeReferenceGradientElement(int i, const VectReal_wp& JacobDof,
						 const VectReal_wp& JacobQuadrature,
						 VectReal_wp& GradJacobQuad) = 0;

  public:
    virtual void ConstructFiniteElement(const string& name_elt) = 0;
    virtual void UpdateInterpolationElement() = 0;
    virtual void ClearFiniteElement() = 0;
    virtual void SplitMeshForParallelComputation(const string& name_elt) = 0;    
    virtual void ComputeNumberOfDofs() = 0;
    virtual void PutOtherGlobalDofs() = 0;
    virtual void CheckContinuity() = 0;

#ifdef MONTJOIE_WITH_TRANSMISSION        
    virtual void PartMeshTransmission() = 0;
    virtual void TreatTransmission(const IVect&) = 0;
    virtual void SendTransmissionDofs(const IVect&, int&, int&, IVect&, IVect&) = 0;
    virtual void DistributeTransmissionDofs(int, int, IVect&, IVect&) = 0;
#endif

    virtual void TreatGibc(const IVect&) = 0;
    virtual void InitGibcReferences(int) = 0;
    
    virtual void
    ComputeEnHnOnBoundary(const MeshInterpolationFEM<Dimension>& interp, const Vector<Real_wp>& U0,
			  VectReal_wp& trace_En, VectReal_wp& trace_Hn,
			  bool assemble = true, bool compute_H = true) const = 0;
    
    virtual void
    ComputeEnHnOnBoundary(const MeshInterpolationFEM<Dimension>& interp, const VectComplex_wp& U0,
			  VectComplex_wp& trace_En, VectComplex_wp& trace_Hn,
			  bool assemble = true, bool compute_H = true) const = 0;    
    
  };

}

#define MONTJOIE_FILE_VAR_GEOMETRY_PROBLEM_HXX
#endif
