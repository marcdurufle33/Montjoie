#ifndef MONTJOIE_FILE_PROD_MAT_VECT_HCURL_3D_HXX

namespace Montjoie
{
  
  template <class T, class T0, class TypeEquation>
  void MltAddHcurl3D(const T& alpha, const GlobalGenericMatrix<T0>& nat_mat,
                     const SeldonTranspose& trans, int level,
                     const FemMatrixFreeClass<T0, TypeEquation>& A,
                     const Vector<T>& B2, const T& beta, Vector<T>& C2, bool assemble);
  
  template <class T, class TypeEquation>
  void MltAddHcurlHdg3D(const T& alpha, const GlobalGenericMatrix<T>& nat_mat,
			const SeldonTranspose& trans, int level,
			const FemMatrixFreeClass<T, TypeEquation>& A,
			const Vector<T>& B2, const T& beta, Vector<T>& C2, bool assemble);
  
  template<class Complexe, class TypeEquation>
  void GetExtrapolationSecondOrder_Sipg(const EllipticProblem<TypeEquation>& var,
                                        int num_elem, const Vector<Complexe>& B,
                                        ExtrapolVariablesProductFEM<Complexe, TypeEquation>& var_extra,
                                        const GlobalGenericMatrix<Complexe>& nat_mat,
                                        const ElementReference<Dimension3, 2>& Fb);

  template<class Complexe, class TypeEquation>
  void GetExtrapolationAndExchange_Sipg(const GlobalGenericMatrix<Complexe>& nat_mat,
                                        const FemMatrixFreeClass<Complexe, TypeEquation>& A,
                                        const Vector<Complexe>& B,
                                        ExtrapolVariablesProductFEM<Complexe, TypeEquation>& var_extra);

  template<class T0, class TypeEquation, class Complexe>
  void MltAdd_ElementHcurlSipg(const T0& alpha, const GlobalGenericMatrix<Complexe>& nat_mat,
                               const FemMatrixFreeClass<Complexe, TypeEquation>& A,
                               int i, const Vector<Complexe>& B, Vector<Complexe>& C,
                               const ElementReference<Dimension3, 2>& Fb);

  template <class T0, class Complexe, class TypeEquation, class Vector1>
  void MltAdd_SipgMaxwell(const T0& alpha, const GlobalGenericMatrix<Complexe>& nat_mat,
                          const SeldonTranspose& trans, int level,
                          const FemMatrixFreeClass<Complexe, TypeEquation>& A,
                          const Vector1& B2, const T0& beta, Vector1& C2, bool assemble);

} // namespace Montjoie

#define MONTJOIE_FILE_PROD_MAT_VECT_HCURL_3D_HXX
#endif
