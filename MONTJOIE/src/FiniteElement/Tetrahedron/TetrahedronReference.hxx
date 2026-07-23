#ifndef MONTJOIE_FILE_TETRAHEDRON_REFERENCE_HXX

#include "TetrahedronGeomReference.hxx"

namespace Montjoie
{
  
  //! base class for finite tetrahedral elements
  template<int type>
  class TetrahedronReference : public ElementReference<Dimension3, type>
  {
  public :
    int type_interpolation; //!< regular or Hesthaven (Lobatto)
    
  protected :
    TetrahedronGeomReference Fb_geom;

  public :
    TetrahedronReference();

    size_t GetMemorySize() const;
    void SetInputData(const string&, const Vector<string>&);

    inline const TetrahedronGeomReference& GetGeometricTetrahedron() const { return Fb_geom; }

    /****************************
     * Initialization functions *
     ****************************/

    
    void ConstructNumberMap(NumberMap& map, int dg) const;
    void ConstructFiniteElement(int r, int rgeom = 0, int rquad = 0,
				int type_quad = -1, int rsurf_tri = 0, int rsurf_quad = 0,
				int type_surf_tri = -1, int type_surf_quad = -1, int gauss_z = -1);
    
  protected :
    void ConstructQuadrature(int r, int type_quad);
            
    template<int t>
    friend ostream& operator <<(ostream& out, const TetrahedronReference<t>& e);
    
  };


} // namespace Montjoie

#define MONTJOIE_FILE_TETRAHEDRON_REFERENCE_HXX
#endif

