#ifndef MONTJOIE_FILE_TETRAHEDRON_HCURL_LOBATTO_HXX

namespace Montjoie
{
  
  //! finite element class for Nedelec's second family on tetrahedra
  class TetrahedronHcurlLobatto : public TetrahedronReference<2>
  {
  protected :
    TetrahedronClassical tet_h1;
    VolumeHcurlReference elt_hcurl;
    
  public :
    TetrahedronHcurlLobatto();

    void ConstructNumberMap(NumberMap& nmap, int dg) const;
    
    // construction of finite element
    void ConstructFiniteElement(int r, int rgeom = 0, int rquad = 0, int type_quad = -1,
				int rsurf_tri = 0, int rsurf_quad = 0,
                                int type_surf_tri = -1, int type_surf_quad = -1, int gauss_z = -1);

    void ComputeValuesPhiRef(const R3&, VectR3&) const;
    void ComputeCurlPhiRef(const R3&, VectR3&) const;

    void GetValuePhiOnQuadraturePoint(int, VectR3& phi) const;
    
    void ApplyCh(const VectReal_wp&, VectReal_wp&) const;
    void ApplyCh(const VectComplex_wp&, VectComplex_wp&) const;
    
  };
  
  ostream& operator <<(ostream& out, const TetrahedronHcurlLobatto& e);

} // namespace Montjoie

#define MONTJOIE_FILE_TETRAHEDRON_HCURL_LOBATTO_HXX
#endif
