#ifndef MONTJOIE_FILE_QUADRANGLE_DG_ORTHO_HXX

namespace Montjoie
{
  
  //! classical nodal finite element on triangles
  /*!
    Basis functions are Lagrange functions generating P_r space
    Interpolation points are Hesthaven's points
   */
  class QuadrangleDgOrtho : public QuadrangleReference<1>
  {
  protected :
    Matrix<Real_wp> LegendrePolynom;
    VectReal_wp InvWeightPolynomial;
    Matrix<int> NumFct2D;
    
  public :
    QuadrangleDgOrtho();

    size_t GetMemorySize() const;
    
    // construction of finite element
    void ConstructFiniteElement(int r, int rgeom = 0, int rquad = 0, int type_quad = -1,
                                int rsurf = 0, int type_surf = -1);

    void ConstructFiniteElement1D(int r, int rgeom, int rquad, int type_quad);

  protected:
    void ConstructFunctions();

    template<class Vector1, class Vector2>
    void ComputeProjectionDofGen(const Vector1& feval, Vector2& contrib) const;
    
  public:
    void ComputeValuesPhiRef(const R2& pointloc, VectReal_wp& phi) const;
    void ComputeGradientPhiRef(const R2& pointloc, VectR2& grad_phi) const;

    void ComputeProjectionDofRef(const VectReal_wp& feval, VectReal_wp& contrib) const;
    void ComputeProjectionDofRef(const VectComplex_wp& feval, VectComplex_wp& contrib) const;
    
  };

}

#define MONTJOIE_FILE_QUADRANGLE_DG_ORTHO_HXX
#endif

