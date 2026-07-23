#ifndef MONTJOIE_FILE_WEDGE_HCURL_FIRST_FAMILY_HXX

namespace Montjoie
{
  
  //! finite element class for Nedelec's first family on triangular prisms
  /*!
    Basis functions are interpolatory functions :
    phi^{hcurl}_i(x, y) phi^GL_j(z) 
    phi^{h1}_i(x, y) phi^G_j(z) e_z
    with phi^G 1-D Lagrange function associated with Gauss points
    phi^GL 1-D Lagrange function associated with Gauss-Lobatto points
    phi^{hcurl} basis functions of edge triangular element
    phi^{h1} basis functions of nodal triangular element
   */
  class WedgeHcurlFirstFamily : public WedgeReference<2>
  {
  protected :
    Matrix<int> CoordinateDofs;
    Globatto<Real_wp> lobz_hcurl, lobz_h1;
    Matrix<Real_wp> GL_GX;
    TriangleGeomReference element_tri_h1;
    
  public :
    WedgeHcurlFirstFamily();
    
    void ConstructNumberMap(NumberMap& map, int dg) const;
    size_t GetMemorySize() const;
    
    // construction of finite element
    void ConstructFiniteElement(int r, int rgeom = 0, int rquad = 0, int type_quad = -1,
				int rsurf_tri = 0, int rsurf_quad = 0,
                                int type_surf_tri = -1, int type_surf_quad = -1, int gauss_z = -1);
    
  protected:
    void ConstructFunctions();

    template<class Vector1, class Vector2>
    void ComputeGaussIntegralSurfaceGen(const Vector1& feval, Vector2& res, int num_loc) const;
    
  public:
    void ComputeValuesPhiRef(const R3& ptloc, VectR3& phi) const;
    void ComputeCurlPhiRef(const R3& ptloc, VectR3& dphi) const;
    
    virtual void ComputeGaussIntegralSurfaceRef(const VectReal_wp & feval,
						VectReal_wp& res, int num_loc) const;
    
    virtual void ComputeGaussIntegralSurfaceRef(const VectComplex_wp & feval,
						VectComplex_wp& res, int num_loc) const;
    
  };
  
  ostream& operator <<(ostream& out, const WedgeHcurlFirstFamily& e);
  
} // namespace Montjoie

#define MONTJOIE_FILE_WEDGE_HCURL_FIRST_FAMILY_HXX
#endif

