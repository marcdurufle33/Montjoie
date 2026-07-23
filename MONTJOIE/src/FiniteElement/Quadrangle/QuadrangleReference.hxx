#ifndef MONTJOIE_FILE_QUADRANGLE_REFERENCE_HXX

#include "QuadrangleGeomReference.hxx"

namespace Montjoie
{
  
  //! Base class for quadrilateral finite element
  /*!
    this class mainly defines the Fi transformation and DFi on quadrilateral elements
   */
  template<int type>
  class QuadrangleReference : public ElementReference<Dimension2, type>
  {
  protected :
    //! quadrature rule for 1-D integration
    Globatto<Real_wp> lob_quad, lob_basis;
    QuadrangleGeomReference Fb_geom;
  
  public :
    Matrix<int> NumQuad2D, CoordinateQuad;
    
    QuadrangleReference();
    
    size_t GetMemorySize() const;

    QuadrangleGeomReference& GetGeometricElement();
    const QuadrangleGeomReference& GetGeometricElement() const;
    
    /****************************
     * Initialization functions *
     ****************************/
        
    // sets map in order to number dofs of the mesh
    void ConstructNumberMap(NumberMap& map, int dg) const;
    void ConstructFiniteElement(int r, int rgeom = 0, int rquad = 0,
                                int type_quad = -1, int rsurf = 0, int type_surf = -1);

    virtual void ConstructFiniteElement1D(int r, int rgeom, int rquad, int type_quad);
    
    void ConstructQuadrature(int r, int type_quad, int rsurf = 0);

    template<int t>
    friend ostream& operator <<(ostream& out,const QuadrangleReference<t>& e);
    
  };

}

#define MONTJOIE_FILE_QUADRANGLE_REFERENCE_HXX
#endif
