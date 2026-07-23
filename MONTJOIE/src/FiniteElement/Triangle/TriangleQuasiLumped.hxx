#ifndef MONTJOIE_FILE_TRIANGLE_QUASI_LUMPED_HXX

namespace Montjoie
{
 
  //! continuous finite elements on triangles with quasi-lumped mass matrix
  class TriangleQuasiLumped : public TriangleClassical
  {
  public :
    TriangleQuasiLumped();
    
    // construction of finite element
    void ConstructFiniteElement(int r, int rgeom = 0, int rquad = 0, int type_quad = -1,
                                int rsurf = 0, int type_surf = -1);
    
  protected:
    template<class Vector1, class Vector2>
    void ApplyChGen(const Vector1& Uh, Vector2& Vh) const;
    
    template<class Vector1, class Vector2>
    void ApplyChTransposeGen(const Vector1& Uh, Vector2& Vh) const;
    
  public:
    void ApplyCh(const VectReal_wp& Vh, VectReal_wp& Uh) const;
    void ApplyCh(const VectComplex_wp& Vh, VectComplex_wp& Uh) const;

    void ApplyChTranspose(const VectReal_wp& Vh, VectReal_wp& Uh) const;
    void ApplyChTranspose(const VectComplex_wp& Vh, VectComplex_wp& Uh) const;

    friend ostream& operator <<(ostream& out, const TriangleQuasiLumped& e);
    
  };

} // namespace Montjoie

#define MONTJOIE_FILE_TRIANGLE_QUASI_LUMPED_HXX
#endif
