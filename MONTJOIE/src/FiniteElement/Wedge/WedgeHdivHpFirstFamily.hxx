#ifndef MONTJOIE_FILE_WEDGE_HDIV_HP_FIRST_FAMILY_HXX

namespace Montjoie
{
  
  //! finite element class for Nedelec's first family on prisms
  /*!
    Hierarchical basis functions are constructed with Jacobi polynomials
    and generate the following finite element space
    D_r(x, y) P_{r-1}(z)  x  P_{r-1}(x, y) P_r(z)    
   */
  class WedgeHdivHpFirstFamily : public WedgeReference<3>
  {
  protected :
    
    Matrix<Real_wp> LegendrePolynom;
    Matrix<Real_wp> JacobiPolynom11;
    Vector<Matrix<Real_wp> > jacobi_2ip1_pol;
    Matrix<Real_wp> ValLeg, ValuePhiQuad2D, ValuePhiTri2D;
    
  public :
    WedgeHdivHpFirstFamily();
    
    void ConstructNumberMap(NumberMap& map, int dg) const;
    size_t GetMemorySize() const;
    
    // construction of finite element
    void ConstructFiniteElement(int r, int rgeom = 0, int rquad = 0, int type_quad = -1,
				int rsurf_tri = 0, int rsurf_quad = 0, int type_surf_tri = -1,
				int type_surf_quad = -1, int gauss_z = -1);
    
  protected:
    void ConstructFunctions();

    template<class Vector1, class Vector2>
    void ComputeProjectionDofGen(const Vector1& feval, Vector2& contrib) const;
 
    template<class Vector1, class Vector2>
    void ComputeGaussIntegralSurfaceGen(const Vector1& feval, Vector2& res, int num_loc) const;
    
  public:
    void ComputeValuesPhiRef(const R3& ptloc, VectR3& phi) const;
    void ComputeDivPhiRef(const R3& ptloc, VectReal_wp& dphi) const;

    virtual void ComputeProjectionDofRef(const VectReal_wp& feval, VectReal_wp& contrib) const;
    virtual void ComputeProjectionDofRef(const VectComplex_wp& feval, VectComplex_wp& contrib) const;
    
    virtual void ComputeGaussIntegralSurfaceRef(const VectReal_wp & feval,
						VectReal_wp& res, int num_loc) const;
    
    virtual void ComputeGaussIntegralSurfaceRef(const VectComplex_wp & feval,
						VectComplex_wp& res, int num_loc) const;

    void ModifySignProjectionSurface(VectReal_wp& contrib, int num_loc) const;
    void ModifySignProjectionSurface(VectComplex_wp& contrib, int num_loc) const;
    
  };
  
  ostream& operator <<(ostream& out, const WedgeHdivHpFirstFamily& e);

} // namespace Montjoie

#define MONTJOIE_FILE_WEDGE_HDIV_HP_FIRST_FAMILY_HXX
#endif

