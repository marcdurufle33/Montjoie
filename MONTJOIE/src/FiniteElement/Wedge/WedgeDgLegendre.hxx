#ifndef MONTJOIE_FILE_WEDGE_DG_LEGENDRE_HXX

namespace Montjoie
{
  
  //! class for orthogonal finite element on wedges with DG formulation
  /*!
    The generated space is P_r(x, y, z)
   */
  class WedgeDgLegendre : public WedgeReference<1>
  {
  protected :
    TinyVector<Matrix<Real_wp, General, RowSparse>, 5> ShLoc;

    Matrix<Real_wp, General, RowSparse> ChX, ChY, ChZ;
    Matrix<Real_wp, General, RowSparse> RhLoc, ShQuad;
    Array3D<int> NumFct3D; 
    VectReal_wp points1d_y, weights1d_y;
    VectReal_wp DerivDxtildeDx, DerivDxtildeDy;
    Matrix<Real_wp> ShDense;
    Matrix<int> NumOrtho2D;
    VectReal_wp InvWeight2D;
    TriangleGeomReference function_tri;
    Globatto<Real_wp> gauss_x, gauss_y;
    
    Matrix<Real_wp> ShTri;
    Matrix<Real_wp, General, RowSparse> DerivShTri, DerivShQuad, DerivShGauss;
    TinyVector<Matrix<Real_wp, General, RowSparse>, 5> DerivShLoc;
    
    Vector<Matrix<Real_wp> > ProjOperatorQuadOrderY;
    
    Vector<Vector<R2> > ProjOperatorDxtildeTri;
    Vector<TinyVector<Vector<R2>, 3> > ProjOperatorDxtildeQuad;
    
    Vector<Matrix<Real_wp, General, RowSparse> > ProjOperatorDerivQuadOrder;
    Vector<Matrix<Real_wp, General, RowSparse> > ProjOperatorDerivTriOrder;
    Vector<Matrix<Real_wp, General, RowSparse> > ProjOperatorDerivQuadOrderY;
    
  public :
    
    WedgeDgLegendre();

    size_t GetMemorySize() const;
    
    // construction of finite element
    void ConstructFiniteElement(int r, int rgeom = 0, int rquad = 0, int type_quad = -1,
				int rsurf_tri = 0, int rsurf_quad = 0,
                                int type_surf_tri = -1, int type_surf_quad = -1, int gauss_z = -1);
    
  protected :
    void ConstructFunctions();
    void ConstructMassMatrix();
    void ConstructStiffnessMatrix();

    template<class Vector1>
    void SolveMassMatrixGen(Vector1& x) const;
    
    template<class Vector1>
    void MltMassMatrixGen(Vector1& y) const;
    
    template<class Vector1, class Vector2>
    void ApplyChGen(const Vector1& U, Vector2& V) const;
    
    template<class Vector1, class Vector2>
    void ApplyChTransposeGen(const Vector1& U, Vector2& V) const;
    
    template<class Vector1, class Vector2>
    void ApplyRhQuadratureGen(const Vector1& U, Vector2& V) const;
    
    template<class Vector1, class Vector2>
    void ApplyRhQuadratureTransposeGen(const Vector1& U, Vector2& V) const;
    
    template<class T0, class Vector1, class Vector2>
    void ApplyShQuadratureGen(const T0& alpha, int num_loc, const Vector1& U,
			      Vector2& V, int r = 0) const;
	  
    template<class Vector1, class Vector2>
    void ApplyShQuadratureTransposeGen(int num_loc, const Vector1& U, Vector2& V, int r = 0) const;
    
    template<class T0, class Vector1, class Vector2>
    void ApplyNablaShQuadratureGen(const T0& alpha, int num_loc, const Vector1& U,
				   Vector2& V, int r = 0) const;
	  
    template<class Vector1, class Vector2>
    void ApplyNablaShQuadratureTransposeGen(int num_loc, const Vector1& U,
					    Vector2& V, int r = 0) const;

    template<class Vector1, class Vector2>
    void ComputeProjectionDofGen(const Vector1& feval, Vector2& contrib) const;

  public :
    void ComputeTriangularInterpolationProjector(const IVect& order_elt,
                                                 const ElementReference<Dimension2, 1>& Fb,
                                                 const Vector<VectR2>& Pts);
    
    void ComputeQuadrangularInterpolationProjector(const IVect& order_elt,
                                                   const ElementReference<Dimension2, 1>& Fb,
                                                   const Vector<VectR2>& Pts);
    
    void ComputeValuesPhiRef(const R3& pointloc, VectReal_wp& ) const;
    void ComputeGradientPhiRef(const R3& pointloc, VectR3& ) const;

    void SolveMassMatrix(VectReal_wp&) const;
    void SolveMassMatrix(VectComplex_wp& x) const;

    void MltMassMatrix(VectReal_wp&) const;
    void MltMassMatrix(VectComplex_wp&) const;

    void ApplyCh(const VectReal_wp& Vh, VectReal_wp& Uh) const;
    void ApplyCh(const VectComplex_wp& Vh, VectComplex_wp& Uh) const;

    void ApplyChTranspose(const VectReal_wp& Vh, VectReal_wp& Uh) const;
    void ApplyChTranspose(const VectComplex_wp& Vh, VectComplex_wp& Uh) const;

    void ApplyRhQuadrature(const VectReal_wp& Uh, VectReal_wp& Vh) const;
    void ApplyRhQuadrature(const VectComplex_wp& Uh, VectComplex_wp& Vh) const;
    
    void ApplyRhQuadratureTranspose(const VectReal_wp& Uh, VectReal_wp& Vh) const;
    void ApplyRhQuadratureTranspose(const VectComplex_wp& Uh, VectComplex_wp& Vh) const;

    void ApplyShQuadratureTranspose(int num_loc, const VectReal_wp& Uh, VectReal_wp& Vh, int r = 0) const;
    void ApplyShQuadratureTranspose(int num_loc, const VectComplex_wp& Uh,
				    VectComplex_wp& Vh, int r = 0) const;
    
    void ApplyShQuadrature(const Real_wp& alpha, int num_loc, const VectReal_wp& Uh,
			   VectReal_wp& Vh, int r = 0) const;

    void ApplyShQuadrature(const Complex_wp& alpha, int num_loc, const VectComplex_wp& Uh,
			   VectComplex_wp& Vh, int r = 0) const;

    void ApplyNablaShQuadratureTranspose(int num_loc, const VectReal_wp& Uh,
					 VectReal_wp& Vh, int r = 0) const;

    void ApplyNablaShQuadratureTranspose(int num_loc, const VectComplex_wp& Uh,
					 VectComplex_wp& Vh, int r = 0) const;
    
    void ApplyNablaShQuadrature(const Real_wp& alpha, int num_loc, const VectReal_wp& Uh,
				VectReal_wp& Vh, int r = 0) const;

    void ApplyNablaShQuadrature(const Complex_wp& alpha, int num_loc, const VectComplex_wp& Uh,
				VectComplex_wp& Vh, int r = 0) const;

    void ComputeProjectionDofRef(const VectReal_wp& feval, VectReal_wp& contrib) const;
    void ComputeProjectionDofRef(const VectComplex_wp& feval, VectComplex_wp& contrib) const;
    
  };

} // namespace Montjoie

#define MONTJOIE_FILE_WEDGE_DG_LEGENDRE_HXX
#endif
