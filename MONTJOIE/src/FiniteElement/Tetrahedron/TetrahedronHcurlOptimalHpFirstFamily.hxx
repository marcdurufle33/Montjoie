#ifndef MONTJOIE_FILE_TETRAHEDRON_HCURL_OPTIMAL_HP_FIRST_FAMILY_HXX

namespace Montjoie
{
  
  //! finite element class for Nedelec's first family on tetrahedra
  /*!
    Hierarchical basis functions are constructed with Jacobi
    polynomials and generate the space Rr (Nedelec's first family space)
   */
  class TetrahedronHcurlOptimalHpFirstFamily : public TetrahedronReference<2>
  {    
  protected :
    
    Matrix<Real_wp> LegendrePolynom;
    Array3D<int> NumDofsX, NumDofsY, NumDofsZ;
    
  public :
    TetrahedronHcurlOptimalHpFirstFamily();

    void ConstructNumberMap(NumberMap& map, int dg) const;
    size_t GetMemorySize() const;
    
    // construction of finite element
    void ConstructFiniteElement(int r, int rgeom = 0, int rquad = 0, int type_quad = -1,
				int rsurf_tri = 0, int rsurf_quad = 0,
                                int type_surf_tri = -1, int type_surf_quad = -1, int gauss_z = -1);

  protected:
    void ConstructFunctions();    
    
    template<class Vector1, class Vector2>
    void ComputeProjectionDofGen(const Vector1& feval, Vector2& contrib) const;

  public:
    void ComputeValuesPhiRef(const R3& ptloc, VectR3& phi) const;
    void ComputeCurlPhiRef(const R3& ptloc, VectR3& dphi) const;
        
    virtual void ComputeProjectionDofRef(const VectReal_wp& feval, VectReal_wp& contrib) const;
    virtual void ComputeProjectionDofRef(const VectComplex_wp& feval, VectComplex_wp& contrib) const;

    void ModifySignProjectionSurface(VectReal_wp& contrib, int num_loc) const;
    void ModifySignProjectionSurface(VectComplex_wp& contrib, int num_loc) const;

  };
  
  ostream& operator <<(ostream& out, const TetrahedronHcurlOptimalHpFirstFamily& e);
  
} // namespace Montjoie

#define MONTJOIE_FILE_TETRAHEDRON_HCURL_OPTIMAL_HP_FIRST_FAMILY_HXX
#endif

