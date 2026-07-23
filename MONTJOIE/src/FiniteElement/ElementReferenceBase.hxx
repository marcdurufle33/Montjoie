#ifndef MONTJOIE_FILE_ELEMENT_REFERENCE_BASE_HXX

namespace Montjoie
{
  //! base class for 2-D or 3-D finite elements
  class ElementReference_Base
  {
  public :
    //! how many informations to display
    int print_level;

    // different proposed numberings
    enum {CONTINUOUS, DISCONTINUOUS, HDG};

    //! Warburton trick to get a mass matrix constant
    static bool use_warburton_trick;
    
  protected:
    //! order of approximation for basis functions (dofs)
    int order;
    //! order of quadrature
    int order_quad;
    
    //! number of quadrature points for volume integrals
    int nb_points_quadrature_inside;
    //! number of quadrature points for boundary integrals
    int nb_points_quadrature_boundaries;    
    //! number of dofs on the edges
    int nb_dof_boundaries;
    //! number of dofs on the element
    int nb_dof_loc;
    //! number of points needed to project a function on finite element space
    int nb_points_dof_inside;
    //! list of dof points on each surface
    Vector<IVect> num_dof_points_surf, num_quad_points_surf;
    
    //! QR factorisation of matrix used to project on degrees of freedom
    Matrix<Real_wp> MatProjectionDof;
    //! QR factorisation of matrix used to project on degrees of freedom
    VectReal_wp tauProjectionDof;
    //! weights associated with dofs
    VectReal_wp weights_dof;

    IVect offset_faceSh; //! offset value for Sh (DG method)
    Vector<int> power_two_face, PosDofOnFace;
    
    Matrix<Real_wp, Symmetric, RowSymPacked> mass_matrix; 
    //!< mass matrix \int \varphi_i \varphi_j
    Matrix<Real_wp, Symmetric, RowSymPacked> mass_matrix_chol;
    
    //! value of phi on quadrature points of the boundary
    TinyVector<Matrix<Real_wp, General, RowSparse>, 6> sparse_const_sh;
    //! value of gradient of phi on quadrature points of the boundary
    TinyVector<Matrix<Real_wp>, 6> const_nabla_sh;
    //! half of integration weight on each face
    TinyVector<VectReal_wp, 6> PoidsFlux;

    //! true for finite elements achieving mass lumping (diagonal or block-diagonal mass matrix)
    bool mass_lumping;
    //! true for finite elements achieving mass lumping on cartesian meshes
    bool mass_lumping_ortho;
    //! true for finite elements providing a diagonal mass matrix
    bool diagonal_mass;
    //! true for finite elements that cannot be used for a continuous Galerkin formulation
    //! but only for a discontinuous one
    bool discontinuous_element;
    //! use Piola transform for basis functions ?
    bool use_piola_transform;
    
    //! booleans used to optimize computations
    bool optimized_mass_matrix, optimized_elem_matrix, sparse_mass_matrix, linear_sparse_mass_matrix,
      use_quadrature_for_sh, use_quadrature_for_rh, use_quadrature_free_sh;

    int type_integration_edge, type_integration_tri, type_integration_quad;

  private:
    const VectReal_wp& weights_quad_inside;
    
  public:
    ElementReference_Base(const VectReal_wp& poids);
    virtual ~ElementReference_Base();
    
    // Inline methods
    bool LumpedMassMatrix() const;
    bool MassLumpingOrthogonalElement() const;
    bool DiagonalMassMatrix() const;
    bool DiscontinuousElement() const;
    bool UsePiolaTransform() const;
    void SetPiolaTransform(bool flag = true);
    bool OptimizedComputationMassMatrix() const;
    bool OptimizedComputationElementaryMatrix() const;
    bool SparseMassMatrix() const;
    bool LinearSparseMassMatrix() const;
    bool UseQuadraturePointsForSh() const;
    bool UseQuadraturePointsForRh() const;
    bool UseQuadratureFreeSh() const;
    int GetMaximalOrderRestriction() const;
        
    int GetOrder() const;    
    int GetQuadratureOrder() const;    
    int GetNbDof() const;
    int GetNbPointsQuadratureInside() const;
    int GetNbPointsQuadBoundaries() const;
    int GetNbPointsUsedForSh() const;
    int GetNbDofBoundaries() const;
    int GetNbPointsDofInside() const;
    int GetNbPointsDofSurface(int num_loc) const;

    int GetPointDofNumber(int num_loc, int j) const;
    int GetQuadNumber(int num_loc, int k) const;
    const Vector<IVect>& GetQuadNumbersBoundary() const;
    int GetOffsetSh(int) const;
    void SetQuadNumbersBoundary(const Vector<IVect>& num_quad);
    void SetDofNumbersBoundary(const Vector<IVect>& num_dof);
    
    int GetTypeIntegrationEdge() const;
    int GetTypeIntegrationTriangle() const;
    int GetTypeIntegrationQuadrangle() const;
    bool IsTangentialDof(int j, int num_loc) const;

    const Real_wp& WeightsDofND(int i) const;
    const VectReal_wp& WeightsDofND() const;
    const VectReal_wp& WeightsND() const;
    const VectReal_wp& GetFluxWeight(int num_loc) const;

    // Inline virtual methods
    //virtual void CancelHighOrderDofs(IVect&, int r);
    
    virtual void ComputeGaussIntegralSurfaceRef(const VectReal_wp & feval,
						VectReal_wp& res, int num_loc) const;
    
    virtual void ComputeGaussIntegralSurfaceRef(const VectComplex_wp & feval,
						VectComplex_wp& res, int num_loc) const;

    virtual void ApplyRh(const VectReal_wp& Uh, VectReal_wp& Vh) const;
    virtual void ApplyRh(const VectComplex_wp& Uh, VectComplex_wp& Vh) const;
    
    virtual void ApplyRhTranspose(const VectReal_wp& Uh, VectReal_wp& Vh) const;
    virtual void ApplyRhTranspose(const VectComplex_wp& Uh, VectComplex_wp& Vh) const;

    virtual void ApplyRhQuadrature(const VectReal_wp& Uh, VectReal_wp& Vh) const;
    virtual void ApplyRhQuadrature(const VectComplex_wp& Uh, VectComplex_wp& Vh) const;
    
    virtual void ApplyRhQuadratureTranspose(const VectReal_wp& Uh, VectReal_wp& Vh) const;
    virtual void ApplyRhQuadratureTranspose(const VectComplex_wp& Uh, VectComplex_wp& Vh) const;
        
    virtual void ApplyRhBoundary(const VectReal_wp& Uh, VectReal_wp& Vh) const;
    virtual void ApplyRhBoundary(const VectComplex_wp& Uh, VectComplex_wp& Vh) const;
    
    virtual void ApplyRhBoundaryTranspose(const VectReal_wp& Uh, VectReal_wp& Vh) const;
    virtual void ApplyRhBoundaryTranspose(const VectComplex_wp& Uh, VectComplex_wp& Vh) const;

    virtual void ApplyShTranspose(int num_loc, const VectReal_wp& Uh, VectReal_wp& Vh, int r = 0) const;
    virtual void ApplyShTranspose(int num_loc, const VectComplex_wp& Uh,
				  VectComplex_wp& Vh, int r = 0) const;
    
    virtual void ApplySh(const Real_wp& alpha, int num_loc, const VectReal_wp& Uh,
			 VectReal_wp& Vh, int r = 0) const;

    virtual void ApplySh(const Complex_wp& alpha, int num_loc, const VectComplex_wp& Uh,
			 VectComplex_wp& Vh, int r = 0) const;
    
    virtual void ApplyShQuadratureTranspose(int num_loc, const VectReal_wp& Uh,
					    VectReal_wp& Vh, int r = 0) const;

    virtual void ApplyShQuadratureTranspose(int num_loc, const VectComplex_wp& Uh,
					    VectComplex_wp& Vh, int r = 0) const;
    
    virtual void ApplyShQuadrature(const Real_wp& alpha, int num_loc, const VectReal_wp& Uh,
				   VectReal_wp& Vh, int r = 0) const;

    virtual void ApplyShQuadrature(const Complex_wp& alpha, int num_loc, const VectComplex_wp& Uh,
 				   VectComplex_wp& Vh, int r = 0) const;

    virtual void ApplyNablaShTranspose(int num_loc, const VectReal_wp& Uh,
				       VectReal_wp& Vh, int r = 0) const;

    virtual void ApplyNablaShTranspose(int num_loc, const VectComplex_wp& Uh,
				       VectComplex_wp& Vh, int r = 0) const;
    
    virtual void ApplyNablaSh(const Real_wp& alpha, int num_loc, const VectReal_wp& Uh,
			      VectReal_wp& Vh, int r = 0) const;

    virtual void ApplyNablaSh(const Complex_wp& alpha, int num_loc, const VectComplex_wp& Uh,
			      VectComplex_wp& Vh, int r = 0) const;
    
    virtual void ApplyNablaShQuadratureTranspose(int num_loc, const VectReal_wp& Uh,
						 VectReal_wp& Vh, int r = 0) const;

    virtual void ApplyNablaShQuadratureTranspose(int num_loc, const VectComplex_wp& Uh,
						 VectComplex_wp& Vh, int r = 0) const;
    
    virtual void ApplyNablaShQuadrature(const Real_wp& alpha, int num_loc, const VectReal_wp& Uh,
					VectReal_wp& Vh, int r = 0) const;

    virtual void ApplyNablaShQuadrature(const Complex_wp& alpha, int num_loc, const VectComplex_wp& Uh,
					VectComplex_wp& Vh, int r = 0) const;

    // other methods
    virtual size_t GetMemorySize() const;
    
    virtual void SolveMassMatrix(VectReal_wp&) const;
    virtual void SolveMassMatrix(VectComplex_wp& x) const;

    virtual void SolveCholesky(const SeldonTranspose& TransA, VectReal_wp&) const;
    virtual void SolveCholesky(const SeldonTranspose& TransA, VectComplex_wp&) const;
    
    virtual void MltMassMatrix(VectReal_wp&) const;
    virtual void MltMassMatrix(VectComplex_wp&) const;
    
    virtual void ComputeMassMatrix(Matrix<Real_wp, Symmetric, RowSymPacked> & A,
				   const VectReal_wp & coef) const;

    virtual void IntegrateMassMatrix(Matrix<Real_wp, Symmetric, RowSymPacked> & A,
				     const VectReal_wp & coef) const;

    virtual void PickNearDofs(int pos, const VectBool& DofUsed, IVect& ListeDof, int nb_dof) const;    

    
    // Purely virtual methods
    virtual void ConstructNumberMap(NumberMap& nmap, int) const = 0;
    
    virtual void ApplyCh(const VectReal_wp& Vh, VectReal_wp& Uh) const = 0;
    virtual void ApplyCh(const VectComplex_wp& Vh, VectComplex_wp& Uh) const = 0;

    virtual void ApplyChTranspose(const VectReal_wp& Vh, VectReal_wp& Uh) const = 0;
    virtual void ApplyChTranspose(const VectComplex_wp& Vh, VectComplex_wp& Uh) const = 0;
        
    virtual void ComputeProjectionDofRef(const VectReal_wp& feval, VectReal_wp& contrib) const = 0;
    virtual void ComputeProjectionDofRef(const VectComplex_wp& feval, VectComplex_wp& contrib) const = 0;

    virtual void ComputeProjectionSurfaceDofRef(const VectReal_wp& feval, VectReal_wp& contrib, int) const = 0;
    virtual void ComputeProjectionSurfaceDofRef(const VectComplex_wp& feval, VectComplex_wp& contrib, int) const = 0;
    
    // template methods    
    template<class Vector1,class Vector2>
    void ComputeIntegralRef(const Vector1 & feval, Vector2& res) const;
    
    template<class Vector1,class Vector2>
    void ComputeIntegralSurfaceRef(const Vector1 & feval,
				    Vector2& res, int num_loc) const;
    
  protected :
    // internal methods
    template<class Vector1,class Vector2>
    void ComputeGaussIntegralSurfaceGen(const Vector1 & feval,
                                        Vector2& res, int num_loc) const;
    
    template<class Vector1, class Vector2>
    void ApplyRhQuadratureGen(const Vector1& Uh, Vector2& Vh) const;
    
    template<class Vector1, class Vector2>
    void ApplyRhQuadratureTransposeGen(const Vector1& Uh, Vector2& Vh) const;
    
    template<class Vector1, class Vector2>
    void ApplyShTransposeGen(int num_loc, const Vector1& Uh, Vector2& Vh, int r = 0) const;
    
    template<class T0, class Vector1, class Vector2>
    void ApplyShGen(const T0& alpha, int num_loc, const Vector1& Uh, Vector2& Vh, int r = 0) const;
    
    template<class Vector1, class Vector2>
    void ApplyShQuadratureTransposeGen(int num_loc, const Vector1& Uh, Vector2& Vh, int r = 0) const;
    
    template<class T0, class Vector1, class Vector2>
    void ApplyShQuadratureGen(const T0& alpha, int num_loc, const Vector1& Uh,
			      Vector2& Vh, int r = 0) const;
    
    template<class Vector1, class Vector2>
    void ApplyNablaShTransposeGen(int num_loc, const Vector1& Uh, Vector2& Vh, int r = 0) const;
    
    template<class T0, class Vector1, class Vector2>
    void ApplyNablaShGen(const T0& alpha, int num_loc, const Vector1& Uh,
			 Vector2& Vh, int r = 0) const;
    
    template<class Vector1, class Vector2>
    void ApplyNablaShQuadratureTransposeGen(int num_loc, const Vector1& Uh,
					    Vector2& Vh, int r = 0) const;
    
    template<class T0, class Vector1, class Vector2>
    void ApplyNablaShQuadratureGen(const T0& alpha, int num_loc, const Vector1& Uh,
				   Vector2& Vh, int r = 0) const;

  };
  
}

#define MONTJOIE_FILE_ELEMENT_REFERENCE_BASE_HXX
#endif

