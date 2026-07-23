#ifndef MONTJOIE_FILE_TRIANGLE_CLASSICAL_HXX

namespace Montjoie
{
  
  //! classical nodal finite element on triangles
  /*!
    Basis functions are Lagrange functions generating P_r space
    Interpolation points are Hesthaven's points
   */
  class TriangleClassical : public TriangleReference<1>
  {
  protected :
    Matrix<Real_wp> GL_G;
    Matrix<Real_wp> InverseBasisVDM;
    VectReal_wp InvWeightFct; Matrix<int> NumFct2D;
    Globatto<Real_wp> basis_phi1D;
    TinyVector<Matrix<R2>, 6> GradientPhi_Boundary;
    
  public :
    TriangleClassical();

    size_t GetMemorySize() const;

    // construction of finite element
    void ConstructFiniteElement(int r, int rgeom = 0, int rquad = 0, int type_quad = -1,
                                int rsurf = 0, int type_surf = -1);
    
    void SetInterpolationPoints(const VectR2& points);
    
    void ComputeValuesPhiRef(const R2& pointloc, VectReal_wp& phi) const;
    void ComputeGradientPhiRef(const R2& pointloc, VectR2& grad_phi) const;
    
  protected :
    void ConstructFunctions();
    
    template<class Vector1, class Vector2>
    void ApplyShTransposeGen(int num_loc, const Vector1& Uh, Vector2& Vh, int r = 0) const;
    
    template<class T0, class Vector1, class Vector2>
    void ApplyShGen(const T0& alpha, int num_loc, const Vector1& Uh, Vector2& Vh, int r = 0) const;
    
    template<class Vector1, class Vector2>
    void ApplyNablaShTransposeGen(int num_loc, const Vector1& Uh, Vector2& Vh, int r = 0) const;
    
    template<class T0, class Vector1, class Vector2>
    void ApplyNablaShGen(const T0& alpha, int num_loc, const Vector1& Uh,
                      Vector2& Vh, int r = 0) const;
    
  public :
    void PickNearDofs(int pos, const VectBool& DofUsed, IVect& ListeDof, int nb_dof) const;

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

    
    friend ostream& operator <<(ostream& out, const TriangleClassical& e);
    
  };
  
  
  //! class for nodal finite element on triangles and with DG formulation
  /*!
    The meshes can be hybrid (triangles + quadrilaterals)
    Quadrature rules on edges are Gauss-Lobatto formulas
  */
  class TriangleLobatto : public TriangleClassical
  {

  protected :
    Matrix<Real_wp> GL_Gquad;
    
  public :
    TriangleLobatto();

    size_t GetMemorySize() const;
    
    // construction of finite element
    void ConstructFiniteElement(int r, int rgeom = 0, int rquad = 0, int type_quad = -1,
                                int rsurf = 0, int type_surf = -1);
    
  protected :
    template<class Vector1, class Vector2>
    void ComputeGaussIntegralSurfaceGen(const Vector1 & feval,
                                        Vector2& res, int num_loc) const;

    
  public :
    void ComputeGaussIntegralSurfaceRef(const VectReal_wp & feval, VectReal_wp& res, int num_loc) const;
    void ComputeGaussIntegralSurfaceRef(const VectComplex_wp & feval, VectComplex_wp& res, int num_loc) const;

  };
  
} // namespace Montjoie

#define MONTJOIE_FILE_TRIANGLE_CLASSICAL_HXX
#endif
