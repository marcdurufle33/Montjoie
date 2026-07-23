#ifndef MONTJOIE_FILE_TETRAHEDRON_QUASI_LUMPED_HXX

namespace Montjoie
{
  
  //! quasi lumped tetrahedra
  class TetrahedronQuasiLumped : public TetrahedronClassical
  {
  public :
    TetrahedronQuasiLumped();
    
    // construction of finite element
    void ConstructFiniteElement(int r, int rgeom = 0, int rquad = 0, int type_quad = -1,
				int rsurf_tri = 0, int rsurf_quad = 0,
                                int type_surf_tri = -1, int type_surf_quad = -1, int gauss_z = -1);
    
    template<class Vector1, class Vector2>
    void ApplyCh(const Vector1& Uh, Vector2& Vh) const;
    
    template<class Vector1, class Vector2>
    void ApplyChTranspose(const Vector1& Uh, Vector2& Vh) const;
    
    friend ostream& operator <<(ostream& out, const TetrahedronQuasiLumped& e);
    
  };

} // namespace Montjoie

#define MONTJOIE_FILE_TETRAHEDRON_QUASI_LUMPED_HXX
#endif
