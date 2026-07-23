#ifndef MONTJOIE_FILE_PYRAMID_HDIV_OTHER_HXX

namespace Montjoie
{
  //! class for facet elements on pyramids
  /*!
    test class, do not use
  */
  class PyramidHdivOther : public PyramidReference<3>
  {
  public :
    
    //! choice of basis functions to use
    int type_basis;
    Vector< TinyVector< MultivariatePolynomial< Real_wp >, 3 > > P;
    Matrix<Real_wp> ValuePhiQuad2D, ValuePhiTri2D;
    //! available choices
    enum {NIGAM_PHILLIPS_1, NIGAM_PHILLIPS_2, GRAGLIA, ZAGLMAYR, GRADINARU, OPTIMAL, OPTIMAL_HP};
    
  public :
    PyramidHdivOther();

    void ConstructOrthogonalBasis(int r);
    
    void ConstructFunctions();
    void ConstructNumberMap(NumberMap& map, int dg) const;
    
    // construction of finite element
    void ConstructFiniteElement(int r, int rgeom = 0, int rquad = 0, int type_quad = -1,
				int rsurf_tri = 0, int rsurf_quad = 0,
                                int type_surf_tri = -1, int type_surf_quad = -1, int gauss_z = -1);
    
    void ComputeValuesPhiRef(const R3& ptloc, VectR3& phi) const;
    void ComputeDivPhiRef(const R3& ptloc, VectReal_wp& dphi) const;
    
    void ComputeValuesPhiOrthoRef(const R3& pointloc, VectR3& phi) const;
    
    template<class Vector1, class Vector2>
    void ComputeGaussIntegralSurfaceGen(const Vector1& feval, Vector2& res, int num_loc) const;
    
    virtual void ComputeGaussIntegralSurfaceRef(const VectReal_wp & feval,
						VectReal_wp& res, int num_loc) const;
    
    virtual void ComputeGaussIntegralSurfaceRef(const VectComplex_wp & feval,
						VectComplex_wp& res, int num_loc) const;

    void ModifySignProjectionSurface(VectReal_wp& contrib, int num_loc) const;
    void ModifySignProjectionSurface(VectComplex_wp& contrib, int num_loc) const;
    
  };
  
  ostream& operator <<(ostream& out, const PyramidHdivOther& e);
  
} // namespace Montjoie

#define MONTJOIE_FILE_PYRAMID_HDIV_OTHER_HXX
#endif

