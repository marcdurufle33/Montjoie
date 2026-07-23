#ifndef MONTJOIE_FILE_TRIANGLE_QUADRATURE_HXX

namespace Montjoie
{

  //! class for integration over the unit triangle  
  class TriangleQuadrature
  {
  public :
    
    enum { QUADRATURE_GAUSS, QUADRATURE_TENSOR, QUADRATURE_MASS_LUMPED, QUADRATURE_QUASI_LUMPED};
    
    static void ConstructQuadrature(int p, VectR2& points, VectReal_wp& weights,
                                    int type_quad = QUADRATURE_GAUSS);
    
    static void GetTensorizedQuadratureRule(int r, VectR2& points2d, VectReal_wp& weights,
                                            int type_quad = QUADRATURE_GAUSS);
    
    static void CompleteTrianglePointsWithSymmetry(VectR2& points2d, VectReal_wp& weights2d,
                                                   int nb_couples, const IVect& NbPoints_couple);
    
  };
  
}

#define MONTJOIE_FILE_TRIANGLE_QUADRATURE_HXX
#endif
