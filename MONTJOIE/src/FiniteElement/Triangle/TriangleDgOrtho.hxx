#ifndef MONTJOIE_FILE_TRIANGLE_DG_ORTHO_HXX

namespace Montjoie
{
  
  //! class for orthogonal finite element on triangles and with DG formulation
  /*!
    The meshes can be hybrid (triangles + quadrilaterals)
    Hierarchical orthogonal functions are used (Dubiner's basis)
    and generate classical polynomial space P_r
  */
  class TriangleDgOrtho : public TriangleReference<1>
  {
  protected :
    Matrix<int> NumFct2D; VectReal_wp InvWeightFct;
    Matrix<Real_wp, General, RowSparse> ChX, ChY, rh_loc;
    VectReal_wp points1d_y, ValGaussX0, ValGaussX1, ValGaussY0;
    VectReal_wp DerivGaussX0, DerivGaussX1, DerivGaussY0;
    Matrix<Real_wp> Gy_G, dGy_G;
    VectReal_wp DerivDxtildeDx, DerivDxtildeDy;
    Globatto<Real_wp> gauss_x, gauss_y;
    Vector<Matrix<Real_wp> > ProjOperatorOrderY;
    int type_quadrature;
    
  public :
    TriangleDgOrtho();
    
    size_t GetMemorySize() const;
    const VectReal_wp& GetInverseWeightFunction() const;
    
    // construction of finite element
    void ConstructFiniteElement(int r, int rgeom = 0, int rquad = 0, int type_quad = -1,
                                int rsurf = 0, int type_surf = -1);
    
    void ConstructFiniteElement1D(int r, int rgeom, int rquad, int type_quad);
    
  protected :
    void ConstructMassMatrix();
    void ConstructStiffnessMatrix();
    
    template<class Vector1, class Vector2>
    void ApplyRhGen(const Vector1& Uh, Vector2& Vh) const;
    
    template<class Vector1, class Vector2>
    void ApplyRhTransposeGen(const Vector1& Uh, Vector2& Vh) const;

    template<class Vector1, class Vector2>
    void ApplyRhQuadratureGen(const Vector1& Uh, Vector2& Vh) const;
    
    template<class Vector1, class Vector2>
    void ApplyRhQuadratureTransposeGen(const Vector1& Uh, Vector2& Vh) const;

    template<class Vector1, class Vector2>
    void ApplyChGen(const Vector1& Uh, Vector2& Vh) const;
    
    template<class Vector1, class Vector2>
    void ApplyChTransposeGen(const Vector1& Uh, Vector2& Vh) const;
    
    template<class Vector1, class Vector2>
    void ApplyShQuadratureTransposeGen(int num_loc, const Vector1& Uh, Vector2& Vh, int r = 0) const;
    
    template<class T0, class Vector1, class Vector2>
    void ApplyShQuadratureGen(const T0& alpha, int num_loc, const Vector1& Uh,
			      Vector2& Vh, int r = 0) const;
    
    template<class Vector1, class Vector2>
    void ApplyNablaShQuadratureTransposeGen(int num_loc, const Vector1& Uh,
					    Vector2& Vh, int r = 0) const;
    
    template<class T0, class Vector1, class Vector2>
    void ApplyNablaShQuadratureGen(const T0& alpha, int num_loc, const Vector1& Uh,
				   Vector2& Vh, int r = 0) const;

    template<class Vector1, class Vector2>
    void ComputeProjectionDofGen(const Vector1& feval, Vector2& contrib) const;

  public :
    void ComputeInterpolationProjectorOrder(const IVect& order_elt,
					    const ElementReference<Dimension1, 1>& Fb,
					    const Vector<VectReal_wp>& Pts);
    
    void ComputeValuesPhiRef(const R2& pointloc, VectReal_wp& ) const;
    void ComputeGradientPhiRef(const R2& pointloc, VectR2& ) const;
    
    void ApplyRh(const VectReal_wp& Vh, VectReal_wp& Uh) const;
    void ApplyRh(const VectComplex_wp& Vh, VectComplex_wp& Uh) const;

    void ApplyRhTranspose(const VectReal_wp& Uh, VectReal_wp& Vh) const;
    void ApplyRhTranspose(const VectComplex_wp& Uh, VectComplex_wp& Vh) const;

    void ApplyRhQuadrature(const VectReal_wp& Uh, VectReal_wp& Vh) const;
    void ApplyRhQuadrature(const VectComplex_wp& Uh, VectComplex_wp& Vh) const;
    
    void ApplyRhQuadratureTranspose(const VectReal_wp& Uh, VectReal_wp& Vh) const;
    void ApplyRhQuadratureTranspose(const VectComplex_wp& Uh, VectComplex_wp& Vh) const;

    void ApplyCh(const VectReal_wp& Vh, VectReal_wp& Uh) const;
    void ApplyCh(const VectComplex_wp& Vh, VectComplex_wp& Uh) const;

    void ApplyChTranspose(const VectReal_wp& Vh, VectReal_wp& Uh) const;
    void ApplyChTranspose(const VectComplex_wp& Vh, VectComplex_wp& Uh) const;
        
    void ComputeProjectionDofRef(const VectReal_wp& feval, VectReal_wp& contrib) const;
    void ComputeProjectionDofRef(const VectComplex_wp& feval, VectComplex_wp& contrib) const;

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
    
  };

}

#define MONTJOIE_FILE_TRIANGLE_DG_ORTHO_HXX
#endif
