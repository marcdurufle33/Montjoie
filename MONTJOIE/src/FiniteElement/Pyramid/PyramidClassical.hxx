#ifndef MONTJOIE_FILE_PYRAMID_CLASSICAL_HXX

namespace Montjoie
{
  
  //! nodal finite element on classic pyramid
  /*!
    Nodal basis functions are based on non-regular interpolation points
    These points are constructed to coincide with Gauss-Lobatto points on the quadrilateral base
    and with Hesthaven's points on triangular faces.
    Nodal functions are constructed by using a Vandermonde matrix
    VDM_{i, j} = psi_i(xi_j)
    where xi_j are interpolation points and psi_i orthogonal functions of the space
    P_r(x, y, z) + \sum_{k=0}^{r-1} P_k(x, y) (x y/(1-z))^(r-k) 
  */
  class PyramidClassical : public PyramidReference<1>
  {
  protected :    
    
    VectReal_wp InvWeightFct;
    Matrix<Real_wp> InverseBasisVDM;
    Matrix<int> NumDofs2D_quad;
    Array3D<int> NumFct3D;
    TriangleClassical* function_basis_tri;
    QuadrangleGauss* function_basis_quad;
    TinyVector<Matrix<R3>, 5> GradientPhi_Boundary;
    Globatto<Real_wp> lob_basis, gauss_basis;
    Vector<Matrix<Real_wp> > ProjOperatorTriSurfOrder;

    VectReal_wp DerivDxtildeDx, DerivDxtildeDz, DerivDytildeDy, DerivDytildeDz;
    Matrix<Real_wp, General, RowSparse> rh_loc;
    Matrix<Real_wp> GL_G;
    int type_quad_basis;
    
  public :
    PyramidClassical();
    size_t GetMemorySize() const;
    
    // construction of finite element
    void ConstructFiniteElement(int r, int rgeom = 0, int rquad = 0, int type_quad = -1,
                                int rtri = 0, int r_quad = 0,
                                int type_surf_tri = -1, int type_surf_quad = -1, int gauss_z = -1);    
  protected :
    void ConstructFunctions();
    void ConstructStiffnessMatrix();

    template<class Vector1, class Vector2>
    void ApplyRhQuadratureGen(const Vector1& Uh, Vector2& Vh) const;
    
    template<class Vector1, class Vector2>
    void ApplyRhQuadratureTransposeGen(const Vector1& Uh, Vector2& Vh) const;
    
    template<class Vector1, class Vector2>
    void ApplyShTransposeGen(int num_loc, const Vector1& Uh, Vector2& Vh, int r = 0) const;
    
    template<class T0, class Vector1, class Vector2>
    void ApplyShGen(const T0& alpha, int num_loc, const Vector1& Uh, Vector2& Vh, int r = 0) const;

    const Matrix<Real_wp>& GetProjectionOperatorFace(int num_loc, int r) const;
    
    template<class Vector1, class Vector2>
    void ApplyNablaShTransposeGen(int num_loc, const Vector1& Uh, Vector2& Vh, int r = 0) const;
    
    template<class T0, class Vector1, class Vector2>
    void ApplyNablaShGen(const T0& alpha, int num_loc, const Vector1& Uh,
			 Vector2& Vh, int r = 0) const;

    template<class Vector1,class Vector2>
    void ComputeGaussIntegralSurfaceGen(const Vector1 & feval, Vector2& res, int num_loc) const;
    
  public :    
    void ComputeValuesPhiRef(const R3& x, VectReal_wp& phi) const;
    void ComputeGradientPhiRef(const R3& x, VectR3& grad_phi) const;

    void ComputeTriangularInterpolationProjector(const IVect& order_elt,
                                                 const ElementReference<Dimension2, 1>&,
                                                 const Vector<VectR2>& Pts);
    
    void ComputeQuadrangularInterpolationProjector(const IVect& order_elt,
                                                   const ElementReference<Dimension2, 1>&,
                                                   const Vector<VectR2>& Pts);
    
    void ApplyRhQuadrature(const VectReal_wp& Uh, VectReal_wp& Vh) const;
    void ApplyRhQuadrature(const VectComplex_wp& Uh, VectComplex_wp& Vh) const;
    
    void ApplyRhQuadratureTranspose(const VectReal_wp& Uh, VectReal_wp& Vh) const;
    void ApplyRhQuadratureTranspose(const VectComplex_wp& Uh, VectComplex_wp& Vh) const;

    void ApplyShTranspose(int num_loc, const VectReal_wp& Uh, VectReal_wp& Vh, int r = 0) const;
    void ApplyShTranspose(int num_loc, const VectComplex_wp& Uh, VectComplex_wp& Vh, int r = 0) const;
    
    void ApplySh(const Real_wp& alpha, int num_loc, const VectReal_wp& Uh, VectReal_wp& Vh, int r = 0) const;
    void ApplySh(const Complex_wp& alpha, int num_loc, const VectComplex_wp& Uh,
		 VectComplex_wp& Vh, int r = 0) const;
    
    void ApplyNablaShTranspose(int num_loc, const VectReal_wp& Uh, VectReal_wp& Vh, int r = 0) const;
    void ApplyNablaShTranspose(int num_loc, const VectComplex_wp& Uh, VectComplex_wp& Vh, int r = 0) const;
    
    void ApplyNablaSh(const Real_wp& alpha, int num_loc, const VectReal_wp& Uh,
		      VectReal_wp& Vh, int r = 0) const;
    
    void ApplyNablaSh(const Complex_wp& alpha, int num_loc, const VectComplex_wp& Uh,
		      VectComplex_wp& Vh, int r = 0) const;

    void ComputeGaussIntegralSurfaceRef(const VectReal_wp & feval, VectReal_wp& res, int num_loc) const;
    void ComputeGaussIntegralSurfaceRef(const VectComplex_wp & feval, VectComplex_wp& res, int num_loc) const;
    
    friend ostream& operator <<(ostream& out,const PyramidClassical& e);
    
  };

  
  //! same class as PyramidClassical except Gauss-Lobatto quadrature on quadrilateral base
  class PyramidLobatto : public PyramidClassical
  {
  public :
    PyramidLobatto();
    
    // construction of finite element
    void ConstructFiniteElement(int r, int rgeom = 0, int rquad = 0, int type_quad = -1,
				int rsurf_tri = 0, int rsurf_quad = 0,
                                int type_surf_tri = -1, int type_surf_quad = -1, int gauss_z = -1);
    
  };

} // namespace Montjoie

#define MONTJOIE_FILE_PYRAMID_CLASSICAL_HXX
#endif
