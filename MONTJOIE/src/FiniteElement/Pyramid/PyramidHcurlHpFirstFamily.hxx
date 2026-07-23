#ifndef MONTJOIE_FILE_PYRAMID_HCURL_HP_FIRST_FAMILY_HXX

namespace Montjoie
{
  
  //! finite element class for Nedelec's first family on pyramid
  /*!
    Hierarchical basis functions are constructed with Jacobi polynomials
    and generate the following finite element space
    C_{r-1}^3 + { x^p y^p (1-z)^p [y, x, xy],  0 <= p <= r-1 }
    + { x^m y^(n+2) (1-z)^(n+1) [1, 0, x],  x^(n+2) y^m (1-z)^(n+1) [0, 1, y], 0 <= m <= n <= r-2 }
    + { x^p y^q (1-z)^r [1, 0, x],  x^q y^p (1-z)^r [0, 1, y],   p <= r-1; q <= r}    
    with C_r = Span( x^i y^j (1-z)^k  i, j <= k <= r )
    Here x, y, z are coordinates on the cube [-1,1]^2 x [0,1] obtained as
    x = xhat / (1-zhat), y = yhat / (1-zhat), z = zhat
    with xhat, yhat, zhat coordinates in the symmetric pyramid
  */
  class PyramidHcurlHpFirstFamily : public PyramidReference<2>
  {
  public:
    int type_basis;
    enum {OPTIMAL, NIGAM_PHILLIPS_1};

  protected :
    
    Matrix<Real_wp> LegendrePolynom;
    Matrix<Real_wp> JacobiPolynom11;
    Vector<Matrix<Real_wp> > JacobiPolynomMij;
    Matrix<Real_wp> OperatorDofTri;
    Array3D<int> NumDofsX, NumDofsY, NumDofsZ;  
    
  public :
    PyramidHcurlHpFirstFamily();
    
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
  
  ostream& operator <<(ostream& out, const PyramidHcurlHpFirstFamily& e);
  
} // namespace Montjoie

#define MONTJOIE_FILE_PYRAMID_HCURL_HP_FIRST_FAMILY_HXX
#endif

