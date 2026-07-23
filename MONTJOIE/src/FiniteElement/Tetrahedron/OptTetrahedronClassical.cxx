#ifndef MONTJOIE_FILE_OPT_TETRAHEDRON_CLASSICAL_CXX

namespace Montjoie
{
    
  template<int r> template<class Vector1, class Vector2>
  inline void OptTetrahedronClassical<r>::ApplyCh(const Vector1& Vh, Vector2& Uh) const
  {
    Mlt(mat_ch, Vh, Uh); 
  }
    
  template<int r> template<class Vector1, class Vector2>
  inline void OptTetrahedronClassical<r>::ApplyChTranspose(const Vector1& Uh, Vector2& Vh) const
  {
    MltTrans(mat_ch, Uh, Vh); 
  }
    
  template<int r> template<class Vector1, class Vector2>
  inline void OptTetrahedronClassical<r>::ApplyRh(const Vector1& Uh, Vector2& Vh) const
  {
    Mlt(mat_rigid, Uh, Vh); 
  }
  
  template<int r> template<class Vector1, class Vector2>
  inline void OptTetrahedronClassical<r>::ApplyRhTranspose(const Vector1& Uh, Vector2& Vh) const
  {
    MltTrans(mat_rigid, Uh, Vh); 
  }

  
  template<int r> template<class Vector1, class Vector2>
  inline void OptTetrahedronClassical<r>::MltMassMatrix(const Vector1& Uh, Vector2& Vh) const
  {
    Mlt(mat_mass, Uh, Vh);
  }
  
}

#define MONTJOIE_FILE_OPT_TETRAHEDRON_CLASSICAL_CXX
#endif
