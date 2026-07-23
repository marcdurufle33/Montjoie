#ifndef MONTJOIE_FILE_QUADRANGLE_QUADRATURE_HXX

namespace Montjoie
{

  //! class for integration over the unit square [0, 1]^2
  class QuadrangleQuadrature
  {
  public :
    
    static void ConstructQuadrature(int p, VectR2& points, VectReal_wp& weights,
                                    int type_quad = Globatto<Real_wp>::QUADRATURE_GAUSS,
				    Real_wp alpha = 0, Real_wp beta = 0);
    
    static void ConstructQuadrilateralNumbering(int order, Matrix<int>& NumNodes2D,
						Matrix<int>& CoordinateNodes);
    
  };
}

#define MONTJOIE_FILE_QUADRANGLE_QUADRATURE_HXX
#endif
