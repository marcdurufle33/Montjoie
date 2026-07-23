#ifndef MONTJOIE_FILE_PROD_MAT_VECT_HCURL_2D_HXX

namespace Montjoie
{
  
  template<class T0, class TypeEquation, class Complexe>
  void MltAdd_ElementHcurl(const T0& alpha, const GlobalGenericMatrix<Complexe>& nat_mat,
			   const FemMatrixFreeClass<Complexe, TypeEquation>& A,
			   int i, const Vector<Complexe>& B, Vector<Complexe>& C,
                           const ElementReference<Dimension2, 2>& Fb);
  
  template <class T0, class Complexe, class TypeEquation>
  void MltAddHcurl2D(const T0& alpha, const GlobalGenericMatrix<Complexe>& nat_mat,
		     const SeldonTranspose& trans, int level,
		     const FemMatrixFreeClass<Complexe, TypeEquation>& A,
		     const Vector<T0>& B2, const T0& beta, Vector<T0>& C2, bool assemble);
  
} // namespace Montjoie

#define MONTJOIE_FILE_PROD_MAT_VECT_HCURL_2D_HXX
#endif
