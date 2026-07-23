#ifndef MONTJOIE_FILE_HEXAHEDRON_REFERENCE_HXX

#include "HexahedronGeomReference.hxx"

namespace Montjoie
{
  
  //! base class for hexahedral finite element
  template<int type>
  class HexahedronReference : public ElementReference<Dimension3, type>
  {
  protected :
    //! Lobatto basis functions
    Globatto<Real_wp> lob_basis, lob_quad;
    
    Array3D<int> NumQuad3D;
    Matrix<int> CoordinateQuad3D, CoordinateQuad2D;
    HexahedronGeomReference Fb_geom;
    
    bool serendip_element;
    
  public :
    HexahedronReference();

    size_t GetMemorySize() const;
    
    void ConstructNumberMap(NumberMap& map, int dg) const;
    void ConstructFiniteElement(int r, int rgeom = 0, int rquad = 0, int type_quad = -1,
                                int rsurf_tri = 0, int rsurf_quad = 0, int type_surf_tri = -1,
                                int type_surf_quad = -1, int gauss_z = -1) = 0;
    
  protected :
    void ConstructQuadrature(int r, int type_quad, int rsurf = 0, int type_quad_surf = -1);
        
    template<int t>
    friend ostream& operator <<(ostream& out, const HexahedronReference<t>& e);
    
  };
  
} // namespace Montjoie

#define MONTJOIE_FILE_HEXAHEDRON_REFERENCE_HXX
#endif

