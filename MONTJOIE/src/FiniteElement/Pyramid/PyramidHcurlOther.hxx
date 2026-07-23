#ifndef MONTJOIE_FILE_PYRAMID_HCURL_OTHER_HXX

namespace Montjoie
{

  // class for Zgainski triangle
  class TriangleHcurlOther : public TriangleReference<2>
  {
  public :
    int type_basis;
    
    void ConstructFiniteElement(int r, int rgeom = 0, int rquad = 0, int type_quad = -1,
				int rsurf = 0, int type_surf = -1);
    
    void ComputeValuesPhiRef(const R2& pt, VectR2& phi) const;
    void ComputeCurlPhiRef(const R2& pt, VectReal_wp& curl_phi) const;

  };


  // class for Zgainski quadrilateral
  class QuadrangleHcurlOther : public QuadrangleReference<2>
  {
  public :
    int type_basis;
    
    void ConstructFiniteElement(int r, int rgeom = 0, int rquad = 0, int type_quad = -1,
				int rsurf = 0, int type_surf = -1);
    
    void ComputeValuesPhiRef(const R2& pt, VectR2& phi) const;
    void ComputeCurlPhiRef(const R2& pt, VectReal_wp& curl_phi) const;
    
  };
    
  //! finite element class for edge elements on pyramids
  /*!
    test class, do not use
   */
  class PyramidHcurlOther : public PyramidReference<2>
  {
  public :
    
    //! choice of basis functions to use
    int type_basis;
    //! available choices
    enum {ZGAINSKI, NIGAM_PHILLIPS_1, NIGAM_PHILLIPS_2, GRAGLIA, ZAGLMAYR, DOUCET, GRADINARU};
    
  public :
    PyramidHcurlOther();

    void ConstructFunctions();
    void ConstructNumberMap(NumberMap& map, int dg) const;
    
    // construction of finite element
    void ConstructFiniteElement(int r, int rgeom = 0, int rquad = 0, int type_quad = -1,
				int rsurf_tri = 0, int rsurf_quad = 0,
                                int type_surf_tri = -1, int type_surf_quad = -1, int gauss_z = -1);
    
    void ComputeValuesPhiRef(const R3& ptloc, VectR3& phi) const;
    void ComputeCurlPhiRef(const R3& ptloc, VectR3& dphi) const;
    
  };
  
  ostream& operator <<(ostream& out, const PyramidHcurlOther& e);
  
} // namespace Montjoie

#define MONTJOIE_FILE_PYRAMID_HCURL_OTHER_HXX
#endif

