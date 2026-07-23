#ifndef MONTJOIE_FILE_WEDGE_REFERENCE_HXX

#include "WedgeGeomReference.hxx"

namespace Montjoie
{
  
  //! base class for finite prismatic elements
  template<int type>
  class WedgeReference : public ElementReference<Dimension3, type>
  {
  protected :
    Matrix<int> NumQuadTri, CoordinateQuad;
    VectReal_wp points1d_z, weights1d_z;
    WedgeGeomReference Fb_geom;
    
  public :
    WedgeReference();

    size_t GetMemorySize() const;    
    
    /****************************
     * Initialization functions *
     ****************************/
    
    void ConstructFiniteElement(int r, int rgeom = 0, int rquad = 0, int type_quad = -1,
				int rsurf_tri = 0, int rsurf_quad = 0,
                                int type_surf_tri = -1, int type_surf_quad = -1, int gauss_z = -1);
    
    void ConstructNumberMap(NumberMap& nmap, int dg) const;

  protected :
    void ConstructQuadrature(int r, int type_quad, int gauss_z,
                             int rsurf_tri = 0, int rsurf_quad = 0,
                             int type_surf_tri = -1, int type_surf_quad = -1);
    
    template<int t>
    friend ostream& operator <<(ostream& out, const WedgeReference<t>& e);
    
  };

} // end namespace

#define MONTJOIE_FILE_WEDGE_REFERENCE_HXX
#endif
