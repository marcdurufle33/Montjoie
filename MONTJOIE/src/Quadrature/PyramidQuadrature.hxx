#ifndef MONTJOIE_FILE_PYRAMID_QUADRATURE_HXX

namespace Montjoie
{

  //! class for integration over the symmetric pyramid  
  class PyramidQuadrature
  {
  public :
    // when the symmetric pyramid is transformed into a cube
    // we have the integral \int (1-z)^2 f(x, y, z) dx dy dz
    // where x, y, z belongs to the cube
    // Because of the presence of (1-z)^2 we can use Gauss-Jacobi rules
    // QUADRATURE_GAUSS : exact formulas on the cube for Q_{2r+1}
    // QUADRATURE_JACOBI1 : exact formulas on the cube for (1-z) Q_{2r+1}
    // QUADRATURE_JACOBI2 : exact formulas on the cube for (1-z)^2 Q_{2r+1}
    // QUADRATURE_NON_PRODUCT : exact formulas on the pyramid for P_{2r}
    enum {QUADRATURE_GAUSS, QUADRATURE_JACOBI1, QUADRATURE_JACOBI2,
          QUADRATURE_NON_PRODUCT};
    
    static void ConstructQuadrature(int r, VectR3& points3d, VectReal_wp& weights3d,
                                    VectReal_wp& points1d_z, VectReal_wp& weights1d_z,
                                    int type_quad);
    
    static void ConstructPolynomialRule(int p, VectR3& points3d, VectReal_wp& weights3d);
    
  };

}

#define MONTJOIE_FILE_PYRAMID_QUADRATURE_HXX
#endif

