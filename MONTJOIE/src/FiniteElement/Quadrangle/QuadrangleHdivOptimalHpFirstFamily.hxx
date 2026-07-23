#ifndef MONTJOIE_FILE_QUADRANGLE_HDIV_OPTIMAL_HP_FIRST_FAMILY_HXX

namespace Montjoie
{
  
  //! class for facet finite element of optimal Nedelec's first family
  /*!
    An hierarchical basis is constructed with Jacobi polynomials
    and generate the optimal first family Q_{r+1, r-1} x Q_{r-1, r+1}
   */
  class QuadrangleHdivOptimalHpFirstFamily : public QuadrangleReference<3>
  {    
  protected :
    Matrix<Real_wp> LegendrePolynom;
    Matrix<Real_wp> JacobiPolynom;
    Matrix<int> NumDofsX, NumDofsY;
    
  public :
    QuadrangleHdivOptimalHpFirstFamily();
    
    // how to number mesh
    void ConstructNumberMap(NumberMap& nmap, int dg) const;
    size_t GetMemorySize() const;
    
    // construction of finite element
    void ConstructFiniteElement(int r, int rgeom = 0, int rquad = 0, int type_quad = -1,
                                int rsurf = 0, int type_surf = -1);

  protected:
    void ConstructFunctions();
    
    template<class Vector1, class Vector2>
    void ComputeProjectionDofGen(const Vector1& feval, Vector2& contrib) const;

  public:
    void ComputeValuesPhiRef(const R2& pointloc, VectR2& ) const;
    void ComputeDivPhiRef(const R2& pointloc, VectReal_wp& ) const;
    
    virtual void ComputeProjectionDofRef(const VectReal_wp& feval, VectReal_wp& contrib) const;  
    virtual void ComputeProjectionDofRef(const VectComplex_wp& feval, VectComplex_wp& contrib) const;

  };
  
  ostream& operator <<(ostream& out, const QuadrangleHdivOptimalHpFirstFamily& e);

} // namespace Montjoie

#define MONTJOIE_FILE_QUADRANGLE_HDIV_OPTIMAL_HP_FIRST_FAMILY_HXX
#endif

