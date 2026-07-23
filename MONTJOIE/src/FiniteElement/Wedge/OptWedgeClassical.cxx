#ifndef MONTJOIE_FILE_OPT_WEDGE_CLASSICAL_CXX

namespace Montjoie
{
  
  template<int r> 
  inline void OptWedgeClassical<r>::ConstructFiniteElement(const WedgeClassical& Fb)
  {
    for (int i = 0; i < Fb.GetNbDof(); i++)
      for (int j = i; j < Fb.GetNbDof(); j++)
	mat_mass(i, j) = Fb.mass_matrix(i, j);
    
    for (int i = 0; i < Fb.GetNbDof(); i++)
      for (int j = 0; j < Fb.GetNbPointsQuadratureInside(); j++)
	{
	  mat_ch(i, j) = Fb.Value_Phi(i, j);
	  for (int k = 0; k < 3; k++)
	    mat_rigid(i, 3*j+k) = Fb.Gradient_Phi(j, i)(k);
	}
  }
    
  
  template<int r> template<class Vector1, class Vector2>
  inline void OptWedgeClassical<r>::ApplyCh(const Vector1& Vh, Vector2& Uh)
  {
    // Mlt(mat_ch, Vh, Uh); 
  }
    
  template<int r> template<class Vector1, class Vector2>
  inline void OptWedgeClassical<r>::ApplyChTranspose(const Vector1& Uh, Vector2& Vh)
  {
    // MltTrans(mat_ch, Vh, Uh); 
  }
    
  template<int r> template<class Vector1, class Vector2>
  inline void OptWedgeClassical<r>::ApplyRh(const Vector1& Uh, Vector2& Vh)
  {
    // Mlt(mat_rigid, Uh, Vh); 
  }
  
  template<int r> template<class Vector1, class Vector2>
  inline void OptWedgeClassical<r>::ApplyRhTranspose(const Vector1& Uh, Vector2& Vh)
  {
    // MltTrans(mat_rigid, Uh, Vh); 
  }

  
  template<int r> template<class Vector1, class Vector2>
  inline void OptWedgeClassical<r>::MltMassMatrix(const Vector1& Uh, Vector2& Vh)
  {
    // Mlt(mat_mass, Uh, Vh);
  }
  
}

#define MONTJOIE_FILE_OPT_WEDGE_CLASSICAL_CXX
#endif
