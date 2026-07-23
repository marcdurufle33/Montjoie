#ifndef MONTJOIE_FILE_VAR_PROBLEM_HXX

namespace Montjoie
{
  
  //! base class with finite element and mesh
  template<class Dimension>
  class VarProblem : public VarComputationProblem,
		     public DistributedProblem<Dimension>,
		     public VarFiniteElement<Dimension>
  {
  public :
    
    // define of R^n, arrays of R^n etc
    typedef typename Dimension::R_N R_N; //!< point in R2 or R3
    typedef typename Dimension::VectR_N VectR_N; //!< vector of points
    typedef typename Dimension::MatrixN_N MatrixN_N; //!< 2x2 or 3x3 matrix
    typedef typename Dimension::DimensionBoundary DimensionBoundary;
    
    //! 2x2 or 3x3 symmetric matrix
    typedef typename Dimension::MatrixN_Nsym MatrixN_Nsym;

  private:
    VarBoundaryCondition_Base& var_boundary;
    
  public :

    template<class TypeEquation>
    VarProblem(EllipticProblem<TypeEquation>&);
            
    // Other methods
    void GetMemoryUsed(map<string, size_t>& var) const;

    void GetInternalNodesElement(int i, int nb_dof_loc, int& nb_dof_edges, int& nb_dof_int,
				 Vector<int>& intern_node) const;

    void GetStaticCondensedRows(IVect& IndexCondensedRows, IVect&, IVect&, IVect&, int& nb_scalar_dof,
				int&, IVect&, Vector<IVect>&) const;
    
    virtual void ConstructFiniteElement(const string& name_elt);
    virtual void UpdateInterpolationElement();
    void ClearFiniteElement();
    
    void CheckContinuity();
    void ComputeDofCoordinates(Vector<R_N>& PointsDof);

  protected:
    void ComputeGeometryQuantity(IVect& OrderFace);

  public:
    void ComputeReferenceGradientElement(int i, const VectReal_wp& JacobDof,
					 const VectReal_wp& JacobQuadrature, VectReal_wp& GradJacobQuad);
    
    void ComputeMassMatrix(bool compute_rho = true, bool delete_points = true);
    int GetMassMatrixType(Vector<bool>&) const;
    int GetElementaryMatrixType() const;
    bool IsSymmetricElementaryMatrix(const GlobalGenericMatrix<Real_wp>&) const;
    bool IsSymmetricElementaryMatrix(const GlobalGenericMatrix<Complex_wp>&) const;
    bool IsSparseElementaryMatrix(const GlobalGenericMatrix<Real_wp>&) const;
    bool IsSparseElementaryMatrix(const GlobalGenericMatrix<Complex_wp>&) const;
    
    void SetInputData(const string& description_field,
		      const VectString& parameters);
    
    void ComputeNumberOfDofs();
    void PutOtherGlobalDofs();

    virtual void
    ComputeEnHnQuadrature(Vector<VectReal_wp>& u_quadrature,
			  Vector<VectReal_wp>& grad_quadrature,
			  int num_elem, const Vector<R_N>& pts, const Vector<R_N>& normale,
			  bool compute_H, VectReal_wp& En_quad, VectReal_wp& Hn_quad) const;

    virtual void
    ComputeEnHnQuadrature(Vector<VectComplex_wp>& u_quadrature,
			  Vector<VectComplex_wp>& grad_quadrature,
			  int num_elem, const Vector<R_N>& pts, const Vector<R_N>& normale,
			  bool compute_H, VectComplex_wp& En_quad, VectComplex_wp& Hn_quad) const;

    void ComputeLocalProlongation(FiniteElementInterpolator& proj,
				  DistributedProblem<Dimension>& var_coarse, int rc, 
				  TinyVector<Matrix<Real_wp>, 4>& ProlongationElement) const;

    void ComputeElementaryMatrix(int, IVect&, VirtualMatrix<Real_wp>&,
				 CondensationBlockSolver_Base<Real_wp>&,
				 const GlobalGenericMatrix<Real_wp>&);
    
    void ComputeElementaryMatrix(int, IVect&, VirtualMatrix<Complex_wp>&,
				 CondensationBlockSolver_Base<Complex_wp>&,
				 const GlobalGenericMatrix<Complex_wp>&);
    
    void CopyInputData(const VarProblem_Base& var);
    
    // Inline methods
    void GetVelocityOnElements(VectReal_wp& velocity);
    Real_wp GetCoefficientPenaltyStiffness(int ref) const;

    void CopyFiniteElement(const DistributedProblem<Dimension>& var);
    
    void SplitMeshForParallelComputation(const string& name_elt);

    int GetNbPointsQuadratureInside(int) const;
    int GetNbDofBoundaries(int num_elem, int nm = 0) const;
    int GetNbLocalDof(int num_elem, int nm = 0) const;
    int GetNbSurfaceDof(int i, int nm = 0) const;
    int GetNbElt() const;
    int GetPrintLevel() const;

    inline VarBoundaryCondition_Base& GetBoundaryConditionProblem() { return var_boundary; }
    inline const VarBoundaryCondition_Base& GetBoundaryConditionProblem() const { return var_boundary; }

    const VectReal_wp& WeightsND(int) const;
    const typename Dimension::DimensionBoundary::VectR_N& PointsQuadratureBoundary(int, int) const;
    const typename Dimension::DimensionBoundary::VectR_N& PointsDofBoundary(int, int) const;
    const ElementGeomReference<Dimension>& GetShapeElement(int i) const;
    const ElementReference_Base& GetReferenceElementBase(int i, int n = 0) const;
    const ElementReference_Base& GetSurfaceElementBase(int i, int n = 0) const;

    virtual bool DiagonalScalarMassMatrix(const ElementReference_Dim<Dimension>&, int i) const;
    virtual bool BlockDiagonalScalarMassMatrix(const ElementReference_Dim<Dimension>&, int i) const;

    VectR_N PointsQuadInsideND(int i) const;
    
    void ComputeEnHnOnBoundary(const MeshInterpolationFEM<Dimension>& interp, const Vector<Real_wp>& U0,
			       VectReal_wp& trace_En, VectReal_wp& trace_Hn,
			       bool assemble = true, bool compute_H = true) const;

    void ComputeEnHnOnBoundary(const MeshInterpolationFEM<Dimension>& interp, const VectComplex_wp& U0,
			       VectComplex_wp& trace_En, VectComplex_wp& trace_Hn,
			       bool assemble = true, bool compute_H = true) const;
    
  };
  
}

#define MONTJOIE_FILE_VAR_PROBLEM_HXX  
#endif
