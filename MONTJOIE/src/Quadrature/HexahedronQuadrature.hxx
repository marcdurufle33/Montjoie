#ifndef MONTJOIE_FILE_HEXAHEDRON_QUADRATURE_HXX

namespace Montjoie
{
  
  //! class for integration over the unit cube [0, 1]^3
  class HexahedronQuadrature
  {
  public :
    
    static void ConstructQuadrature(int p, VectR3& points, VectReal_wp& weights,
                                    int type_quad = Globatto<Real_wp>::QUADRATURE_GAUSS,
				    Real_wp alpha = 0, Real_wp beta = 0);
    
    static void GetEdgeCube(Matrix<int>& hex_edge);
    
    static void ConstructHexahedralNumbering(int order, Array3D<int>& NumNodes3D,
					     Matrix<int>& CoordinateNodes);
    
  };
}

#define MONTJOIE_FILE_HEXAHEDRON_QUADRATURE_HXX
#endif
