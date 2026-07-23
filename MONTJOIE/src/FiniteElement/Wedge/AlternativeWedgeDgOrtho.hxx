#ifndef MONTJOIE_FILE_WEDGE_DG_ORTHO_HXX

namespace Montjoie
{
  
  //! class for nodal finite element on pyramids and with DG formulation
  class WedgeDgOrtho : public WedgeReference<1>
  {
    friend class FiniteElementH1<Dimension3>;
    
  protected :
    Matrix<Real_wp, General, RowSparse> ShLoc, ShLocTranspose;
    TinyVector<Matrix<Real_wp>, 2> ShSurf, ShSurfTranspose;
    Matrix<Real_wp, General, RowSparse> ChX, ChTransposeX, ChY, ChTransposeY, ChZ, ChTransposeZ;
    Matrix<Real_wp, General, RowSparse> RhLoc, RhLocTranspose, ShQuad, ShQuadTranspose;
    Array3D<int> NumQuad3D; 
    VectReal_wp points1d_y, weights1d_y;
    VectReal_wp DerivDxtildeDx, DerivDxtildeDy;
    Matrix<Real_wp> ShDense, ShDenseTranspose;
    Matrix<Real_wp> alpha, beta, gamma;
    Array3D<Real_wp> C, CmO,CpO;
        
  public :
    //! volumetric projector
    typedef DenseProjector<Dimension3> ProjectionOperator;
    //! surfacic projector
    typedef DenseProjector<Dimension2> ProjectionSurfacicOperator;
    
    static const int ORDER_MAX_RESTRICTION_STD = 0; //!< matrix is stored
    static const int ORDER_MAX_MLT_STD = 0; //!< matrix is stored
    
    typedef Dimension3 Dimension; //!< dimension 3
    
    WedgeDgOrtho();
    
    // construction of finite element
    void ConstructQuadrature(int r);
    void ConstructFiniteElement(int r, int rgeom = 0, int rquad = 0, int type_quad = -1,
				int rsurf_tri = 0, int rsurf_quad = 0,
                                int type_surf_tri = -1, int type_surf_quad = -1, int gauss_z = -1);
    
    void ConstructFunctions();
    void ConstructMassMatrix();
    void ConstructStiffnessMatrix();
	  
    void ComputeMassMatrix(Matrix<Real_wp, Symmetric, SymColSkyLine> & A,
                           const VectReal_wp & coef) const;
    void ComputeMassMatrix(Matrix<Real_wp, Symmetric, ArrayRowSymSparse> & A,
                           const VectReal_wp & coef) const;
    void ComputeMassMatrix(Matrix<Real_wp, Symmetric, RowSymPacked> & A,
                           const VectReal_wp & coef) const;
    
    template<class Matrix1>
    void PerformComputationMassMatrix(Matrix1& A, const VectReal_wp& coef) const;
    
    template<class Vector1>
    void SolveMassMatrix(Vector1& x) const;
    
    template<class Vector1>
    void MltMassMatrix(Vector1& y) const;
    
    Real_wp GetMassMatrix(int i, int j) const;
    
    template<class Vector1, class Vector2>
    void ApplyCh(const Vector1& U, Vector2& V) const;
    
    template<class Vector1, class Vector2>
    void ApplyChTranspose(const Vector1& U, Vector2& V) const;
    
    template<class Vector1, class Vector2>
    void ApplyRhQuadrature(Vector1& U, Vector2& V) const;
    
    template<class Vector1, class Vector2>
    void ApplyRhQuadratureTranspose(const Vector1& U, Vector2& V) const;
    
    template<class T0, class Vector1, class Vector2>
    void ApplyShQuadrature(const T0& alpha, const Vector1& U, Vector2& V) const;
	  
    template<class Vector1, class Vector2>
    void ApplyShQuadratureTranspose(const Vector1& U, Vector2& V) const;
    
    void ComputeValuesPhiRef(const R3& pointloc, VectReal_wp& ) const;
    void ComputeGradientPhiRef(const R3& pointloc, VectR3& ) const;
    
    Real_wp GetValuePhi(int num_dof, int num_point) const;
    Real_wp GetValuePhiOnBoundary(int k, int j, int) const;
    
    const R3& GetGradientPhi(int num_dof, int num_point) const;
    const R3& GetGradientPhiOnBoundary(int node, int num_point, int num_loc) const;
    
    template<class Vector1, class Vector2>
    void ComputeProjectionDofRef(const Vector1& feval, Vector2& contrib) const;
    
    template<class Vector1,class Vector2>
    void ComputeIntegralRef(const Vector1 & feval, Vector2 & res) const;
    
    template<class Vector1,class Vector2>
    void ComputeIntegralGradientRef(const Vector1 & feval, Vector2 & res) const;
    
    template<class Vector1,class Vector2>
    void ComputeGaussIntegralSurfaceRef(const Vector1 & feval,
                                        Vector2 & res, int num_loc) const;
    
    template<class Vector1,class Vector2>
    void ComputeIntegralSurfaceRef(const Vector1 & feval,
				   Vector2 & res, int num_loc) const;
    
    template<class Vector1,class Vector2>
    void ComputeIntegralSurfaceGradientRef(const Vector1 & feval,
					   Vector2 & res, int num_loc) const;
    
  };

} // namespace Montjoie

#define MONTJOIE_FILE_WEDGE_DG_ORTHO_HXX
#endif
