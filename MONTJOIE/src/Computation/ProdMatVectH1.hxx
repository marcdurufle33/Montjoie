#ifndef MONTJOIE_FILE_PROD_MAT_VECT_H1_HXX

namespace Montjoie
{

  // functions in ProdMatVectH1.cxx
  template<class TypeEquation, class Vector1, class Vector2, int p, class Vector4, class Vector5>
  void GetExtrapolationU(const EllipticProblem<TypeEquation>& var, int num_elem,
			 const Vector1& Uloc, TinyVector<Vector2, p>& Uquad, Vector4& extrapolU,
                         const ElementReference<typename TypeEquation::Dimension, 1>& Fb,
			 Vector5& Uneighbor, bool compute_uquad);
  
  template<class TypeEquation, class Complexe, class Vector1>
  void GetExtrapolationSecondOrder(const EllipticProblem<TypeEquation>& var,
                                   int num_elem, const Vector1& B,
				   ExtrapolVariablesProductFEM<Complexe, TypeEquation>& var_extra,
                                   const GlobalGenericMatrix<Complexe>& nat_mat,
				   const ElementReference<typename TypeEquation::Dimension, 1>& Fb);
  
  template<class Complexe, class TypeEquation, class T1>
  void GetExtrapolationAndExchange(const GlobalGenericMatrix<Complexe>& nat_mat,
                                   const FemMatrixFreeClass<Complexe, TypeEquation>& A,
                                   const Vector<T1>& B,
				   ExtrapolVariablesProductFEM<T1, TypeEquation>& var_extra);

  template<class Complexe, class TypeEquation, class T0, class Vector1>
  void AddLocalFluxNeighbor(const T0& alpha, const GlobalGenericMatrix<Complexe>& nat_mat,
                            const FemMatrixFreeClass<Complexe, TypeEquation>& A,
                            int i, const Vector1& B, Vector1& C,
                            ExtrapolVariablesProductFEM<Complexe, TypeEquation>& var_extra,
			    const ElementReference<typename TypeEquation::Dimension, 1>& Fb);

  template<class Complexe, class TypeEquation, class T0, class Vector1>
  void CompleteExchangeAndFlux(const T0& alpha, const GlobalGenericMatrix<Complexe>& nat_mat,
                               const FemMatrixFreeClass<Complexe, TypeEquation>& A,
                               const Vector1& B, Vector1& C,
                               ExtrapolVariablesProductFEM<Complexe, TypeEquation>& var_extra);

  template<class T0, class TypeEquation, class Complexe, class Vector1>
  void MltAdd_ElementH1(const T0& alpha, const GlobalGenericMatrix<Complexe>& nat_mat,
			const FemMatrixFreeClass<Complexe, TypeEquation>& A,
			int i, const Vector1& B, Vector1& C,
                        ExtrapolVariablesProductFEM<Complexe, TypeEquation>& var_extra,
			const ElementReference<typename TypeEquation::Dimension, 1>& Fb);

  template <class T0, class Complexe, class Vector1, class TypeEquation>
  void MltAddVectorH1(const T0& alpha, const GlobalGenericMatrix<Complexe>& nat_mat,
		      const SeldonTranspose& trans, int level,
		      const FemMatrixFreeClass<Complexe, TypeEquation>& A,
		      const Vector1& B2, const T0& beta, Vector1& C2, bool assemble);

  template<class Complexe, class Vector2>
  void ExtractApplyDF_PML(const Vector<Complexe>&, int nb_unknowns, int j,
                          bool ortho_z, Dimension2& dim, int offset_elt,
                          const Matrix2_2& dfjm1, const Complexe& coef,
                          Vector2& Vh, Vector2& Vtilde, Vector2& Vterce);

  template<class Complexe, class Vector2>
  void ExtractApplyDF_PML(const Vector<Complexe>& B, int nb_unknowns, int j,
                          bool ortho_z, Dimension3& dim, int offset_elt,
                          const Matrix3_3& dfjm1, const Complexe& coef,
                          Vector2& Vh, Vector2& Vtilde, Vector2& Vterce);

  template<class T0, class Complexe, class TypeEquation,
	   class Vector, class Vector2>
  void MltAdd_ElementH1_Scalar(const T0& alpha,
                               const FemMatrixFreeClass<Complexe, TypeEquation>& A,
			       int i, const Vector& B, Vector& C, Vector2& Uh,
			       Vector2& Vh, const ElementReference<typename TypeEquation::Dimension, 1>& Fb);

  template<class T0, class Complexe, class TypeEquation,
	   class Vector, class Vector2, class Vector3>
  void MltAdd_ElementH1_Vectorial(const T0& alpha,
                                  const FemMatrixFreeClass<Complexe, TypeEquation>& A,
				  int i, const Vector& B, Vector& C, Vector2& dU_loc, Vector3& Uh,
				  Vector3& Vh, const ElementReference<typename TypeEquation::Dimension, 1>& Fb);

  template <class T0, class Vector1, class Complexe, class TypeEquation>
  void MltAdd_SquareHex_VectorialH1(const T0& alpha, int level,
				    const FemMatrixFreeClass<Complexe, TypeEquation>& A,
				    const Vector1& B, const T0& beta, Vector1& C);

  template <class T0, class Vector1, class Complexe, class TypeEquation>
  void MltAdd_SquareHex_ScalarH1(const T0& alpha, int level,
				 const FemMatrixFreeClass<Complexe, TypeEquation>& A,
				 const Vector1& B, const T0& beta, Vector1& C);
  
  // Functions in ProdMatVectScalarH1.cxx  
  template<class Complexe, class TypeEquation, class T1>
  void GetExtrapolationSecondOrder_Scalar(const EllipticProblem<TypeEquation>& var,
                                          int num_elem, const Vector<T1>& B,
					  ExtrapolVariablesProductFEM<T1, TypeEquation>& var_extra,
                                          const GlobalGenericMatrix<Complexe>& nat_mat,
					  const ElementReference<typename TypeEquation::Dimension, 1>& Fb);

  template<class Complexe, class TypeEquation, class T1>
  void GetExtrapolationAndExchange_Scalar(const GlobalGenericMatrix<Complexe>& nat_mat,
                                          const FemMatrixFreeClass<Complexe, TypeEquation>& A,
                                          const Vector<T1>& B,
					  ExtrapolVariablesProductFEM<T1, TypeEquation>& var_extra);

  template<class Complexe, class TypeEquation, class T0, class T1>
  void CompleteExchangeAndFlux_Scalar(const T0& alpha, const GlobalGenericMatrix<Complexe>& nat_mat,
                                      const FemMatrixFreeClass<Complexe, TypeEquation>& A,
                                      Vector<T1>& C, ExtrapolVariablesProductFEM<T1,
                                      TypeEquation>& var_extra);

  template<class T0, class TypeEquation, class Complexe, class T1>
  void MltAdd_ElementH1_Scalar(const T0& alpha, const GlobalGenericMatrix<Complexe>& nat_mat,
			       const FemMatrixFreeClass<Complexe, TypeEquation>& A,
			       int i, const Vector<T1>& B, Vector<T1>& C,
                               ExtrapolVariablesProductFEM<T1, TypeEquation>& var_extra,
			       const ElementReference<typename TypeEquation::Dimension, 1>& Fb);

  template <class T0, class Complexe, class T1, class TypeEquation>
  void MltAddScalarH1(const T0& alpha, const GlobalGenericMatrix<Complexe>& nat_mat,
		      const SeldonTranspose& trans, int level,
		      const FemMatrixFreeClass<Complexe, TypeEquation>& A,
		      const Vector<T1>& B2, const T0& beta, Vector<T1>& C2, bool assemble);

} // namespace Montjoie

#define MONTJOIE_FILE_PROD_MAT_VECT_H1_HXX
#endif
