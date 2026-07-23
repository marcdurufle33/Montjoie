#ifndef MONTJOIE_FILE_TETRAHEDRON_QUADRATURE_HXX

namespace Montjoie
{

  //! class for integration over the unit tetrahedron  
  class TetrahedronQuadrature
  {
  public :
    
    enum { QUADRATURE_GAUSS, QUADRATURE_TENSOR, QUADRATURE_MASS_LUMPED, QUADRATURE_QUASI_LUMPED};
    
    static void ConstructQuadrature(int p, VectR3& points, VectReal_wp& weights,
                                    int type_quad = QUADRATURE_GAUSS);
    
    static void GetTensorizedQuadratureRule(int r, VectR3& points3d, VectReal_wp& weights,
                                            int type_quad = QUADRATURE_GAUSS);
    
    static void CompleteTetrahedronPointsWithSymmetry(VectR3& points3d, VectReal_wp& weights3d,
                                                      int nb_couples,
                                                      const IVect& NbPoints_couple);
    
  };
}

#define MONTJOIE_FILE_TETRAHEDRON_QUADRATURE_HXX
#endif
