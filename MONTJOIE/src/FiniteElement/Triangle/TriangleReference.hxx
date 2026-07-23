#ifndef MONTJOIE_FILE_TRIANGLE_REFERENCE_HXX

#include "TriangleGeomReference.hxx"

namespace Montjoie
{
  
  //! base class for finite elements on triangles
  /*
    this class mainly defines the Fi transformation and DFi
  */
  template<int type>
  class TriangleReference : public ElementReference<Dimension2, type>
  {
  public :
    //! regular interpolation or Lobatto-like interpolation
    static int type_interpolation;
    
  protected:
    TriangleGeomReference Fb_geom;
    
  public :
    
    TriangleReference();

    size_t GetMemorySize() const;    

    TriangleGeomReference& GetGeometricElement();
    const TriangleGeomReference& GetGeometricElement() const;
    
    /****************************
     * Initialization functions *
     ****************************/
    
    
    void ConstructNumberMap(NumberMap& nmap, int) const;    
    void ConstructFiniteElement(int r, int rgeom = 0, int rquad = 0, int type_quad = -1,
                                int rsurf = 0, int type_surf = -1);
    
    void ConstructQuadrature(int r, int type_quad = 0, int rsurf = 0, int type_surf = -1);

    virtual void ConstructFiniteElement1D(int r, int rgeom, int rquad, int type_quad);
    
    template<int t>
    friend ostream& operator <<(ostream& out, const TriangleReference<t>& e);
    
  };
  
  template<int type>
  int TriangleReference<type>::type_interpolation;
  
} // namespace Montjoie

#define MONTJOIE_FILE_TRIANGLE_REFERENCE_HXX
#endif
