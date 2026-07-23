#ifndef MONTJOIE_FILE_ELEMENTARY_MATRIX_HCURL_2D_HXX

namespace Montjoie
{

  template<class Vect1, class Vect2, class Matrix1, class Complexe, class TypeEquation>
  void ComputeVariableElementaryMatrix(int iquad, const Vect1& B, const Vect2& A, Matrix1& mat,
				       const ElementReference<Dimension2, 2>& Fb,
                                       const GlobalGenericMatrix<Complexe>& nat_mat,
                                       const EllipticProblem<TypeEquation>& vars);

  template<class Complexe, class TypeEquation>
  void ComputeElementaryMatrix(int iquad, IVect& num_dof,
                               VirtualMatrix<Complexe>& mat_interac,
			       const GlobalGenericMatrix<Complexe>& nat_mat,
			       const EllipticProblem<TypeEquation>& vars,
                               const ElementReference<Dimension2, 2>& Fb);
  
} // end namespace

#define MONTJOIE_FILE_ELEMENTARY_MATRIX_HCURL_2D_HXX
#endif
