#ifndef MONTJOIE_FILE_WEDGE_HDIV_OPTIMAL_FIRST_FAMILY_HXX

namespace Montjoie
{
  
  //! finite element class for optimal Nedelec's first family on prisms
  /*!
    Basis functions are interpolatory functions :
    phi^{hdiv}_i(x, y) phi^GLint_j(z) 
    phi^{L2}_i(x, y) phi^GL_j(z) e_z
    with phi^GLint 1-D Lagrange function associated with interior Gauss-Lobatto points
    phi^GL 1-D Lagrange function associated with Gauss-Lobatto points
    phi^{hdiv} basis functions of facet triangular element
    phi^{L2} basis functions of nodal triangular element
   */
  class WedgeHdivOptimalFirstFamily : public WedgeReference<3>
  {
  protected :
    VectReal_wp InvWeightBasisRr;
    TriangleHdivFirstFamily Fb_tri;
    Matrix<Real_wp> InverseBasisVDM, LegendrePolynom;
    TriangleGeomReference Fb_triH1;
    Matrix<int> NumDofsXY, NumDofsZ;
    Globatto<Real_wp> lob_lobatto, lob_gauss;
      
  public :
    WedgeHdivOptimalFirstFamily();

    void ConstructNumberMap(NumberMap& map, int dg) const;
    size_t GetMemorySize() const;
    
    // construction of finite element
    void ConstructOrthogonalBasis(int r, int rgeom = 0, int rquad = 0, int type_quad = -1);
    void ConstructFiniteElement(int r, int rgeom = 0, int rquad = 0, int type_quad = -1,
				int rsurf_tri = 0, int rsurf_quad = 0,
                                int type_surf_tri = -1, int type_surf_quad = -1, int gauss_z = -1);

    void ComputeValuesPhiOrthoRef(const R3& ptloc, VectR3& phi) const;
    void ComputeDivPhiOrthoRef(const R3& ptloc, VectReal_wp& dphi) const;
        
  protected:
    void ConstructFunctions();
    
    template<class Vector1, class Vector2>
    void ComputeGaussIntegralSurfaceGen(const Vector1& feval, Vector2& res, int num_loc) const;
    
  public:
    void ComputeValuesPhiRef(const R3& ptloc, VectR3& phi) const;
    void ComputeDivPhiRef(const R3& ptloc, VectReal_wp& dphi) const;
    
    virtual void ComputeGaussIntegralSurfaceRef(const VectReal_wp & feval,
						VectReal_wp& res, int num_loc) const;
    
    virtual void ComputeGaussIntegralSurfaceRef(const VectComplex_wp & feval,
						VectComplex_wp& res, int num_loc) const;
    
    void ModifySignProjectionSurface(VectReal_wp& contrib, int num_loc) const;
    void ModifySignProjectionSurface(VectComplex_wp& contrib, int num_loc) const;
    
  };
  
  ostream& operator <<(ostream& out, const WedgeHdivOptimalFirstFamily& e);

} // namespace Montjoie

#define MONTJOIE_FILE_WEDGE_HDIV_OPTIMAL_FIRST_FAMILY_HXX
#endif

