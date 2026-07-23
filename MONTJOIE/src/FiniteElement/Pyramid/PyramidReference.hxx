#ifndef MONTJOIE_FILE_PYRAMID_REFERENCE_HXX

#include "PyramidGeomReference.hxx"

namespace Montjoie
{
  
  //! base class for finite pyramidal elements
  template<int type>
  class PyramidReference : public ElementReference<Dimension3, type>
  {
  public :    
    int type_interpolation; //!< Hesthaven (Lobatto)
    
  protected :
    VectReal_wp points1d_z, weights1d_z;
    PyramidGeomReference Fb_geom;
    
  public :
    PyramidReference();

    size_t GetMemorySize() const;

    /****************************
     * Initialization functions *
     ****************************/
    
    void ConstructNumberMap(NumberMap& nmap, int dg) const;
    void ConstructFiniteElement(int r, int rgeom = 0, int rquad = 0,
                                int type_quad = -1, int rsurf_tri = 0, int rsurf_quad = 0,
                                int type_surf_tri = -1, int type_surf_quad = -1, int gauss_z = -1);
    
  protected:
    void ConstructQuadrature(int r, int type_quad, int rsurf_tri = 0, int rsurf_quad = 0,
                             int type_surf_tri = -1, int type_surf_quad = -1);
    
    template<int t>
    friend ostream& operator <<(ostream& out, const PyramidReference<t>& e);
    
  };
  
} // namespace Montjoie

#define MONTJOIE_FILE_PYRAMID_REFERENCE_HXX
#endif

